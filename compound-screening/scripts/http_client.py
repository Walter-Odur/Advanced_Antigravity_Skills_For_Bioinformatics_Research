"""Unified HTTP client with rate limiting, retries, and backoff.

Provides `HttpClient` — a single entry-point that combines:

* **Per-API rate limiting** via `_RateLimiter` (cross-process file-lock).
* **Automatic retries** on transient errors (HTTP 429, 5xx).
* **Exponential backoff** with optional *jitter*.
* **Retry-After header** support (server-directed backoff).
* **X-Throttling-Control** proactive backpressure (PubChem / NCBI).
* **Configurable timeouts** per request.

Transport is `urllib.request` (stdlib) — no third-party dependencies.

Typical usage:

```
import http_client

# Scoped to a specific base URL
api_client = http_client.HttpClient(
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
  qps=3
)

# Simple GET with relative path (returns parsed JSON).
data = api_client.fetch_json("esummary.fcgi?db=pubmed&id=123456")

# POST with a JSON body.
data = api_client.fetch_json(
  "esearch.fcgi",
  method="POST",
  json_body={"db": "pubmed", "term": "cancer"},
)

# Raw bytes (e.g. file / PDF download).
pdf = api_client.fetch_bytes(
  "efetch.fcgi?db=pubmed&id=123456&rettype=abstract",
  timeout=60
)
```
"""

import contextlib
import datetime
import email.utils
import gzip
import http.client
import json
import logging
import os
import random
import re
import tempfile
import time
from typing import Any, Iterator
import urllib.error
import urllib.parse
import urllib.request
import urllib.robotparser

# Cross-platform file locking
try:
    import fcntl
    def _lock_file(f):
        fcntl.flock(f, fcntl.LOCK_EX)
    def _unlock_file(f):
        fcntl.flock(f, fcntl.LOCK_UN)
except ImportError:
    import msvcrt
    def _lock_file(f):
        msvcrt.locking(f.fileno(), msvcrt.LK_LOCK, 1)
    def _unlock_file(f):
        try:
            msvcrt.locking(f.fileno(), msvcrt.LK_UNLCK, 1)
        except OSError:
            pass

__all__ = [
    "HttpClient",
    "HttpError",
    "HttpResponse",
    "DEFAULT_BACKOFF_BASE_SECS",
    "DEFAULT_BACKOFF_MAX_SECS",
    "DEFAULT_JITTER_SECS",
    "DEFAULT_MAX_RETRIES",
    "DEFAULT_TIMEOUT_SECS",
    "DEFAULT_USER_AGENT",
    "RETRYABLE_STATUS_CODES",
]

RETRYABLE_STATUS_CODES: frozenset[int] = frozenset({429, 500, 502, 503, 504})
DEFAULT_TIMEOUT_SECS: float = 60.0
DEFAULT_MAX_RETRIES: int = 7
DEFAULT_BACKOFF_BASE_SECS: float = 3.0
DEFAULT_BACKOFF_MAX_SECS: float = 180.0
DEFAULT_JITTER_SECS: float = 0.5
DEFAULT_USER_AGENT: str = os.environ.get(
    "SCIENCE_SKILLS_USER_AGENT",
    "",
)
DEFAULT_CHARSET: str = "utf-8"
PROJECT_NAME: str = "science-skills"

# Regex for parsing X-Throttling-Control status entries.
# Matches e.g. "Request Count status: Red (82%)".
_THROTTLE_STATUS_RE = re.compile(
    r"(\w[\w ]*?) status:\s*(Green|Yellow|Red|Black)\s*\((\d+)%\)"
)

# Backpressure delays (seconds) keyed by throttle colour.
_THROTTLE_BACKPRESSURE: dict[str, float] = {
    "Green": 0.0,
    "Yellow": 1.0,
    "Red": 5.0,
    "Black": 30.0,
}


class _RateLimiter:
  """Enforces a minimum interval between requests.

  Uses a shared lock file so multiple processes respect the same rate limit (set
  lock_file to a path, e.g. /tmp/alphaskills-ncbi.nlm.nih.gov.lock).

  Example:
    limiter = _RateLimiter('ncbi.nlm.nih.gov', qps=10)
  """

  def __init__(self, hostname: str, qps: float):
    """Initialize rate limiter.

    Args:
      hostname: Hostname of the API.
      qps: Maximum queries per second.
    """
    self.hostname = hostname
    self.min_interval = 1.0 / qps
    self.lock_file = os.path.join(tempfile.gettempdir(), f"{PROJECT_NAME}-{hostname}.lock")

  def wait(self, min_sleep: float = 0.0):
    """Block until the next request is allowed.

    Args:
      min_sleep: Additional minimum sleep in seconds.  The actual delay is
        `max(rate_limit_gap, min_sleep)`, which lets callers fold retry back-off
        into the same call.
    """
    with open(self.lock_file, "a+") as f:
      _lock_file(f)
      try:
        f.seek(0)
        content = f.read().strip()
        last_ts = float(content) if content else 0.0
        now = time.monotonic()
        gap = self.min_interval - (now - last_ts)
        delay = max(gap, min_sleep)
        if delay > 0:
          time.sleep(delay)
        f.seek(0)
        f.truncate()
        f.write(str(time.monotonic()))
        f.flush()
      finally:
        _unlock_file(f)


class _RobotsChecker:
  """Lazy-loading robots.txt compliance checker with file-based caching.

  Fetches and caches the robots.txt for a domain on first use.  The parsed
  content is persisted to a per-hostname file under ``/tmp/`` with a 24-hour
  TTL so that short-lived skill processes (the common case) don't re-fetch on
  every invocation.

  If the file cannot be fetched (404, network error, timeout), all URLs are
  allowed — this is the standard behaviour per RFC 9309 §2.4.
  """

  _CACHE_TTL_SECS: float = 86400.0  # 24 hours.

  def __init__(
      self,
      base_url: str,
      user_agent: str,
      timeout: float = 10.0,
      cache_dir: str = None,
  ):
    self._user_agent = user_agent
    self._timeout = timeout
    if cache_dir is None:
      cache_dir = tempfile.gettempdir()
    parsed = urllib.parse.urlparse(base_url)
    self._robots_url = f"{parsed.scheme}://{parsed.netloc}/robots.txt"
    self._cache_file = f"{cache_dir}/{PROJECT_NAME}-{parsed.hostname}.robots"
    self._parser: urllib.robotparser.RobotFileParser | None = None
    self._loaded = False

  def _load(self):
    """Load robots.txt from cache or network (idempotent)."""
    if self._loaded:
      return
    self._loaded = True

    # Try the on-disk cache first.
    cached = self._read_cache()
    if cached is not None:
      parser = urllib.robotparser.RobotFileParser()
      parser.parse(cached.splitlines())
      self._parser = parser
      logging.debug("Loaded robots.txt from cache %s", self._cache_file)
      return

    # Cache miss or stale — fetch from network.
    try:
      req = urllib.request.Request(
          self._robots_url,
          headers={"User-Agent": self._user_agent},
      )
      with urllib.request.urlopen(req, timeout=self._timeout) as resp:
        body = resp.read().decode("utf-8", errors="replace")
      parser = urllib.robotparser.RobotFileParser()
      parser.parse(body.splitlines())
      self._parser = parser
      self._write_cache(body)
      logging.debug("Fetched robots.txt from %s", self._robots_url)
    except Exception:
      # Fail open: missing/broken robots.txt means everything is allowed.
      # Write an empty cache so we don't retry for another 24h.
      self._write_cache("")
      logging.debug(
          "Could not fetch %s — allowing all requests.", self._robots_url
      )

  def _read_cache(self) -> str | None:
    """Read cached robots.txt if present and fresh (< TTL)."""
    try:
      stat = os.stat(self._cache_file)
      age = time.time() - stat.st_mtime
      if age > self._CACHE_TTL_SECS:
        return None
      with open(self._cache_file, "r") as f:
        return f.read()
    except OSError:
      return None

  def _write_cache(self, content: str) -> None:
    """Write robots.txt content to the cache file."""
    try:
      with open(self._cache_file, "w") as f:
        f.write(content)
    except OSError:
      logging.debug("Could not write robots cache to %s", self._cache_file)

  def is_allowed(self, url: str) -> bool:
    """Return True if *url* is allowed by robots.txt for our user-agent."""
    self._load()
    if self._parser is None:
      return True  # Fail open.
    return self._parser.can_fetch(self._user_agent, url)


def _maybe_decompress(
    response: http.client.HTTPResponse,
) -> http.client.HTTPResponse | gzip.GzipFile:
  """Wrap the response in a gzip decompressor if Content-Encoding says so.

  Args:
    response: The response to wrap. Must support .read() and .headers.

  Returns:
    A `GzipFile` wrapper if the response is gzip-encoded, otherwise the
    original response unchanged.
  """
  encoding = response.headers.get("Content-Encoding", "").lower()
  if encoding in ("gzip", "x-gzip"):
    return gzip.GzipFile(fileobj=response)
  return response


class HttpError(Exception):
  """Raised when an HTTP request fails with a non-retryable or exhausted error.

  Attributes:
    status_code: HTTP status code (`None` for network-level failures).
    body: Raw error response body bytes, if available.
    url: The URL that was requested.
  """

  def __init__(
      self,
      message: str,
      *,
      status_code: int | None = None,
      body: bytes | None = None,
      url: str | None = None,
  ):
    super().__init__(message)
    self.status_code = status_code
    self.body = body
    self.url = url

  def json(self) -> Any:
    """Attempt to parse the error body as JSON."""
    # Let the JSON decode fail with JSONDecodeError if the body is empty
    body = self.body or b""
    return json.loads(body.decode("utf-8"))


class HttpResponse:
  """Lightweight wrapper around an HTTP response.

  Attributes:
    data: Raw response body bytes (decompressed if the server sent
      `Content-Encoding: gzip`).
    status_code: HTTP status code.
    headers: Response headers as a dict.
    url: The final URL after any redirects.
    encoding: Character encoding parsed from the `Content-Type` header.  Falls
      back to `"utf-8"` when absent.
  """

  __slots__ = ("data", "status_code", "headers", "url", "encoding")

  def __init__(
      self,
      data: bytes,
      status_code: int,
      headers: dict[str, str],
      url: str,
      encoding: str | None = None,
  ):
    self.data = data
    self.status_code = status_code
    self.headers = headers
    self.url = url
    self.encoding = encoding or "utf-8"

  def json(self) -> Any:
    """Parse the response body as JSON."""
    return json.loads(self.data.decode(self.encoding))

  @property
  def text(self) -> str:
    """Decode the response body using the detected charset."""
    return self.data.decode(self.encoding)

  def __repr__(self) -> str:
    return (
        f"HttpResponse(status={self.status_code}, "
        f"url={self.url!r}, size={len(self.data)})"
    )


def _parse_retry_after(headers) -> float | None:
  """Extract `Retry-After` value from HTTP response headers.

  Handles both formats defined in RFC 7231 §7.1.3:

  * **delta-seconds** — `Retry-After: 120`
  * **HTTP-date** — `Retry-After: Mon, 31 Mar 2026 15:10:00 GMT`

  Args:
    headers: Response headers (dict-like or `http.client.HTTPMessage`).

  Returns:
    Delay in seconds, or `None` if the header is absent or not parseable.
  """
  value = headers.get("Retry-After")
  if value is None:
    return None

  try:
    return float(value)
  except ValueError:
    pass

  try:
    retry_dt = email.utils.parsedate_to_datetime(value)
    delta = (
        retry_dt - datetime.datetime.now(datetime.timezone.utc)
    ).total_seconds()
    return max(0.0, delta)
  except (TypeError, ValueError, OverflowError):
    logging.warning("Failed to parse Retry-After header value: %r", value)
    return None


def _parse_throttle_control(headers) -> float:
  """Parse `X-Throttling-Control` header for proactive backpressure.

  The header (used by PubChem / NCBI) reports real-time usage across three
  dimensions:

      X-Throttling-Control: Request Count status: Green (12%),
          Request Time status: Yellow (55%), Service status: Green (8%)

  Each dimension has a colour: Green (<50%), Yellow (50–75%), Red (>75%), Black
  (blocked).  We return the worst-case backpressure delay across all dimensions
  so the client can slow down **before** hitting a hard 429/503.

  Args:
    headers: Mapping of response headers.

  Returns:
    Recommended additional delay in seconds (0.0 if Green or header
    is absent).
  """
  value = headers.get("X-Throttling-Control")
  if not value:
    return 0.0

  max_delay = 0.0
  for match in _THROTTLE_STATUS_RE.finditer(value):
    colour = match.group(2)
    delay = _THROTTLE_BACKPRESSURE.get(colour, 0.0)
    if delay > max_delay:
      max_delay = delay
  return max_delay


class HttpClient:
  """Rate-limited HTTP client with automatic retries and backoff.

  Uses `urllib.request` as the transport layer. Handles gzip decompression,
  charset detection, and streaming iteration internally.

  Example:

      chembl_api = HttpClient(
          "https://www.ebi.ac.uk/chembl/api/data/",
          qps=5,
          jitter=0.5,
      )
      data = chembl_api.fetch_json("molecule/CHEMBL25.json")
  """

  def __init__(
      self,
      base_url: str,
      qps: float,
      *,
      default_headers: dict[str, str] | None = None,
      max_retries: int = DEFAULT_MAX_RETRIES,
      timeout: float = DEFAULT_TIMEOUT_SECS,
      backoff_base: float = DEFAULT_BACKOFF_BASE_SECS,
      backoff_max: float = DEFAULT_BACKOFF_MAX_SECS,
      jitter: float = DEFAULT_JITTER_SECS,
      user_agent: str = DEFAULT_USER_AGENT,
      retryable_status_codes: frozenset[int] = RETRYABLE_STATUS_CODES,
      permit_direct_access: bool = False,
      robots_cache_dir: str = "/tmp",
  ):
    """Rate-limited HTTP client with automatic retries and backoff.

    Args:
      base_url: Base URL of the API (e.g. "https://eutils.ncbi.nlm.nih.gov/").
      qps: Maximum queries per second (steady-state).
      default_headers: Default HTTP headers to include in every request.
      max_retries: Maximum retry attempts for transient errors.  The total
        number of attempts is `max_retries + 1`.
      timeout: Default per-request timeout in seconds.
      backoff_base: Base delay (seconds) for exponential backoff.
      backoff_max: Cap on backoff delay (seconds).
      jitter: Maximum random jitter (seconds) added uniformly to each backoff
        delay.
      user_agent: Default `User-Agent` header value.  Falls back to the
        `SCIENCE_SKILLS_USER_AGENT` env var.
      retryable_status_codes: HTTP status codes that trigger a retry.
      permit_direct_access: If False (default), fetches and validates requests
        comply with target domain's robots.txt.  Requests to disallowed paths
        raise `HttpError`. If the robots.txt cannot be fetched, all URLs are
        allowed (fail-open). When True, we allow direct access without
        validating against robots.txt.
      robots_cache_dir: Directory for caching robots.txt files (internal testing
        only).
    """
    self.base_url = base_url
    parsed = urllib.parse.urlparse(base_url)
    if not parsed.scheme or not parsed.netloc:
      raise ValueError(f"base_url must be an absolute URL: {base_url!r}")
    self.hostname = parsed.hostname
    self.max_retries = max_retries
    self.timeout = timeout
    self.backoff_base = backoff_base
    self.backoff_max = backoff_max
    self.jitter = jitter
    self.user_agent = user_agent
    self.retryable_status_codes = retryable_status_codes
    self.default_headers = default_headers or {}
    self._limiter = _RateLimiter(self.hostname, qps=qps)
    self._next_min_sleep = 0.0
    self._robots = (
        _RobotsChecker(base_url, self.user_agent, cache_dir=robots_cache_dir)
        if not permit_direct_access
        else None
    )

  def wait(self, min_sleep: float = 0.0) -> None:
    """Wait for the rate limiter without making a request.

    Useful for non-`fetch` operations that still need to respect the
    cross-process rate limit — for example, polling loops or pre-streaming
    handshakes.

    Args:
      min_sleep: Minimum time to sleep in seconds before returning.
    """
    self._limiter.wait(min_sleep=min_sleep)

  def _compute_backoff(
      self, attempt: int, retry_after: float | None = None
  ) -> float:
    """Compute the delay before the next retry.

    Uses exponential backoff `base * 2^attempt` capped at `backoff_max`, with
    optional uniform jitter.  If the server provided a `Retry-After` value, the
    returned delay is at least that large.

    Args:
      attempt: Zero-based retry attempt number (0 = first *retry*).
      retry_after: Optional server `Retry-After` value in seconds.

    Returns:
      Delay in seconds.
    """
    delay = self.backoff_base * (2**attempt)
    if retry_after is not None:
      delay = max(delay, retry_after)
    delay = min(delay, self.backoff_max)
    if self.jitter > 0:
      delay += random.uniform(0, self.jitter)
    return delay

  def _resolve_url(self, url: str) -> str:
    """Resolve *url* against `base_url`."""
    if "://" not in url:
      return urllib.parse.urljoin(self.base_url, url)
    if not url.startswith(self.base_url):
      raise ValueError(f"URL {url!r} does not match base_url {self.base_url!r}")
    return url

  def _build_request(
      self,
      url: str,
      method: str,
      headers: dict[str, str] | None,
      data: bytes | None,
      json_body: Any | None,
  ) -> urllib.request.Request:
    """Build a `urllib.request.Request` from the given parameters."""
    merged_headers = {
        "User-Agent": self.user_agent,
        "Accept-Encoding": "gzip",
    }
    merged_headers.update(self.default_headers)
    if headers:
      merged_headers.update(headers)

    body: bytes | None = data
    if json_body is not None:
      body = json.dumps(json_body).encode("utf-8")
      merged_headers.setdefault("Content-Type", "application/json")

    req = urllib.request.Request(
        url,
        data=body,
        headers=merged_headers,
        method=method,
    )
    return req

  def fetch(
      self,
      url: str,
      *,
      method: str = "GET",
      headers: dict[str, str] | None = None,
      data: bytes | None = None,
      json_body: Any | None = None,
      timeout: float | None = None,
      max_retries: int | None = None,
  ) -> HttpResponse:
    """Execute an HTTP request with rate limiting and retries.

    On each attempt the client:

    1. Waits for the rate limiter (cross-process file-lock).
    2. Sends the request via `urllib.request`.
    3. On success (2xx), returns an `HttpResponse`.
    4. On a retryable HTTP error (429, 5xx) or a network error, sleeps for an
       exponential backoff delay (with optional jitter and `Retry-After`
       support) before retrying.
    5. On a non-retryable HTTP error, raises `HttpError` immediately.

    Args:
      url: Request URL.
      method: HTTP method (`GET`, `POST`, etc.).
      headers: Extra HTTP headers (merged with the default User-Agent).
      data: Raw request body bytes (mutually exclusive with *json_body*).
      json_body: JSON-serializable request body.  Automatically sets
        `Content-Type: application/json`.
      timeout: Per-request timeout in seconds (overrides the client-level
        default).
      max_retries: Override for the maximum number of retry attempts.

    Returns:
      `HttpResponse` containing the response data.

    Raises:
      HttpError: On non-retryable HTTP errors or after exhausting all
          retry attempts.
      ValueError: If both *data* and *json_body* are provided.
    """
    with self._open_stream(
        url, method, headers, data, json_body, timeout, max_retries
    ) as resp:
      stream = _maybe_decompress(resp)
      body = stream.read()
      encoding = resp.headers.get_content_charset() or DEFAULT_CHARSET
      return HttpResponse(
          data=body,
          status_code=resp.status,
          headers=dict(resp.headers),
          url=resp.url,
          encoding=encoding,
      )

  def fetch_json(self, url: str, **kwargs) -> Any:
    """Fetch a URL and parse the response as JSON.

    Convenience wrapper around `fetch()` that adds an
    `Accept: application/json` header (if not already set) and returns the
    parsed JSON body.

    Args:
      url: URL to fetch.
      **kwargs: Keyword arguments to pass to `fetch()`.

    Returns:
      Parsed JSON (dict, list, str, etc.).

    Raises:
      HttpError: On HTTP or network errors.
      json.JSONDecodeError: If the response body is not valid JSON.
    """
    hdrs = kwargs.pop("headers", None) or {}
    hdrs.setdefault("Accept", "application/json")
    resp = self.fetch(url, headers=hdrs, **kwargs)
    return resp.json()

  def fetch_bytes(self, url: str, **kwargs) -> bytes:
    """Fetch a URL and return the raw response body.

    Convenience wrapper for binary downloads (PDFs, images, archives, etc.).

    Args:
      url: URL to fetch.
      **kwargs: Keyword arguments to pass to `fetch()`.

    Returns:
      Raw response bytes.

    Raises:
      HttpError: On HTTP or network errors.
    """
    return self.fetch(url, **kwargs).data

  def fetch_text(self, url: str, **kwargs) -> str:
    """Fetch a URL and return the response body as a decoded string.

    Convenience wrapper for text-based APIs (XML, TSV, plain text).

    Args:
      url: URL to fetch.
      **kwargs: Keyword arguments to pass to `fetch()`.

    Returns:
      Response body decoded using the charset from Content-Type.

    Raises:
      HttpError: On HTTP or network errors.
    """
    return self.fetch(url, **kwargs).text

  @contextlib.contextmanager
  def _open_stream(
      self,
      url: str,
      method: str,
      headers: dict[str, str] | None,
      data: bytes | None,
      json_body: Any | None,
      timeout: float | None,
      max_retries: int | None = None,
  ) -> Iterator[http.client.HTTPResponse]:
    """Open an HTTP response with rate limiting and retries (internal).

    Handles argument validation, rate limiting, request dispatch, and
    error checking.  Yields the raw `http.client.HTTPResponse` and
    guarantees `response.close()` on exit.

    The **connection phase** (before any data flows) is retried on
    transient errors (429, 5xx) with the same backoff logic as
    `fetch()`.  Once a 2xx response is yielded, no further retries
    are attempted — streaming data is not idempotently resumable.

    On **2xx responses** that include an `X-Throttling-Control`
    header, proactive backpressure is applied so the next request
    is delayed before hitting a hard limit.

    Args:
      url: Request URL.
      method: HTTP method.
      headers: Extra HTTP headers.
      data: Raw request body bytes.
      json_body: JSON-serializable request body.
      timeout: Per-request timeout override.
      max_retries: Override for maximum connection retry attempts.

    Yields:
      An open `http.client.HTTPResponse` ready for streaming reads.

    Raises:
      HttpError: On HTTP errors (non-2xx status) or network errors.
      ValueError: If both *data* and *json_body* are provided.
    """
    if data is not None and json_body is not None:
      raise ValueError("Cannot specify both 'data' and 'json_body'.")

    url = self._resolve_url(url)

    # robots.txt compliance check.
    if self._robots and not self._robots.is_allowed(url):
      raise HttpError(
          f"Blocked by robots.txt: {url}",
          url=url,
      )

    effective_timeout = timeout if timeout is not None else self.timeout
    effective_retries = (
        max_retries if max_retries is not None else self.max_retries
    )

    last_exc: Exception | None = None
    next_min_sleep = 0.0
    for attempt in range(effective_retries + 1):
      current_min_sleep = max(next_min_sleep, self._next_min_sleep)
      self._limiter.wait(min_sleep=current_min_sleep)
      self._next_min_sleep = 0.0
      req = self._build_request(url, method, headers, data, json_body)

      try:
        response = urllib.request.urlopen(req, timeout=effective_timeout)
      except urllib.error.HTTPError as exc:
        status = exc.code
        error_body = exc.read()
        retry_after = _parse_retry_after(exc.headers)
        exc.close()

        if (
            status in self.retryable_status_codes
            and attempt < effective_retries
        ):
          next_min_sleep = self._compute_backoff(attempt, retry_after)
          logging.info(
              "HttpClient[%s]: HTTP %d from %s — retrying in ≥%.1fs"
              " (attempt %d/%d)",
              self.hostname,
              status,
              url,
              next_min_sleep,
              attempt + 1,
              effective_retries + 1,
          )
          last_exc = HttpError(
              f"HTTP Error {status} while fetching {url}",
              status_code=status,
              body=error_body,
              url=url,
          )
          continue

        hint = ""
        if status == 403:
          hint = (
              f" (Hint: this may be caused by the User-Agent"
              f" '{self.user_agent}'. Override it by setting the environment"
              f" variable SCIENCE_SKILLS_USER_AGENT, e.g.:"
              f' SCIENCE_SKILLS_USER_AGENT="<enter_your_custom_user_agent>"'
              f" python3 script.py ...)"
          )
        raise HttpError(
            f"HTTP Error {status} while fetching {url}{hint}",
            status_code=status,
            body=error_body,
            url=url,
        ) from exc
      except (urllib.error.URLError, OSError) as exc:
        if attempt < effective_retries:
          next_min_sleep = self._compute_backoff(attempt)
          logging.info(
              "HttpClient[%s]: Network error (%s) — retrying in ≥%.1fs"
              " (attempt %d/%d)",
              self.hostname,
              exc,
              next_min_sleep,
              attempt + 1,
              effective_retries + 1,
          )
          last_exc = exc
          continue
        raise HttpError(
            f"Network error fetching {url}: {exc}", url=url
        ) from exc

      # 2xx success.
      throttle_delay = _parse_throttle_control(response.headers)
      if throttle_delay > 0:
        logging.info(
            "HttpClient[%s]: X-Throttling-Control backpressure %.1fs from %s",
            self.hostname,
            throttle_delay,
            url,
        )
        self._next_min_sleep = throttle_delay

      try:
        yield response
      finally:
        response.close()
      return

    # Should not be reachable.
    raise HttpError(
        f"Max retries ({effective_retries}) exceeded for {url}",
        url=url,
    ) from last_exc

  def stream_lines(
      self,
      url: str,
      *,
      method: str = "GET",
      headers: dict[str, str] | None = None,
      data: bytes | None = None,
      json_body: Any | None = None,
      timeout: float | None = None,
      max_retries: int | None = None,
  ) -> Iterator[str]:
    """Stream an HTTP response line-by-line without buffering.

    Useful for large result sets (e.g. UniProt `/stream` which can
    return up to 10M lines).  The response is streamed via
    `urllib.request` with automatic gzip decompression.

    Rate-limits before each attempt.  The **connection phase** is
    retried on transient errors (429, 5xx), but once data starts
    streaming, no further retries are attempted.

    Args:
      url: Request URL.
      method: HTTP method.
      headers: Extra HTTP headers.
      data: Raw request body bytes (mutually exclusive with *json_body*).
      json_body: JSON-serializable request body.
      timeout: Per-request timeout in seconds (overrides the client default).
        This is the timeout for the initial connection, not for reading
        individual lines.
      max_retries: Override for maximum connection retry attempts.

    Yields:
      Each non-empty line of the response body, decoded as text.

    Raises:
      HttpError: On HTTP errors (non-2xx status).
    """
    with self._open_stream(
        url, method, headers, data, json_body, timeout, max_retries
    ) as response:
      stream = _maybe_decompress(response)
      encoding = response.headers.get_content_charset() or DEFAULT_CHARSET
      for raw_line in stream:
        line = raw_line.decode(encoding).rstrip("\r\n")
        yield line

  def stream_bytes(
      self,
      url: str,
      *,
      method: str = "GET",
      headers: dict[str, str] | None = None,
      data: bytes | None = None,
      json_body: Any | None = None,
      timeout: float | None = None,
      chunk_size: int = 8192,
      max_retries: int | None = None,
  ) -> Iterator[bytes]:
    """Stream an HTTP response as raw byte chunks without buffering.

    Symmetric with `stream_lines` but for binary content (PDFs, archives,
    images, etc.).  Each iteration yields a chunk of up to *chunk_size* bytes.
    Chunks are **transfer-decoded**: if the server applied
    `Content-Encoding: gzip` in transit, the yielded bytes are the
    decompressed content.

    Rate-limits before each attempt.  The **connection phase** is retried on
    transient errors (429, 5xx), but once data starts streaming, no further
    retries are attempted.

    Example:

        with open("paper.pdf", "wb") as f:
            for chunk in client.stream_bytes(url):
                f.write(chunk)

    Args:
      url: Request URL.
      method: HTTP method.
      headers: Extra HTTP headers.
      data: Raw request body bytes (mutually exclusive with *json_body*).
      json_body: JSON-serializable request body.
      timeout: Per-request timeout in seconds (overrides the client default).
        This is the timeout for the initial connection, not for reading
        individual chunks.
      chunk_size: Maximum number of bytes per yielded chunk.
      max_retries: Override for maximum connection retry attempts.

    Yields:
      Non-empty `bytes` chunks of the response body.

    Raises:
      HttpError: On HTTP errors (non-2xx status).
    """
    with self._open_stream(
        url, method, headers, data, json_body, timeout, max_retries
    ) as response:
      stream = _maybe_decompress(response)
      while True:
        chunk = stream.read(chunk_size)
        if not chunk:
          break
        yield chunk
