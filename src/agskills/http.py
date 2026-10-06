"""Rate-limited HTTP client with retries, backoff and a hard deadline.

Transport is ``urllib.request`` from the standard library, so the core of
this package installs without ``requests``.

Differences from the version this replaces, each fixing an observed
failure:

* **Cross-platform cache directory.** The old client cached ``robots.txt``
  under a hardcoded ``/tmp``, which does not exist on Windows. Every write
  failed silently, so ``robots.txt`` was re-fetched on *every* client
  construction - an extra network round trip per command.
* **A real User-Agent.** The old default was the empty string whenever
  ``SCIENCE_SKILLS_USER_AGENT`` was unset. NCBI and several EBI endpoints
  throttle or reject blank agents, producing confusing 403s.
* **Wall-clock rate limiting.** The cross-process limiter persisted
  ``time.monotonic()`` to a shared file. Monotonic clocks have no common
  epoch across processes or reboots, so the stored timestamp was not
  comparable to a fresh reading. It now stores ``time.time()``.
* **The lock is not held while sleeping.** The old limiter slept inside the
  file lock, serialising every concurrent process behind the slowest one.
* **A bounded total deadline.** Old defaults allowed 7 retries with backoff
  capped at 180 s - over six minutes of blocking for a single dead host,
  which is what made the original test suite hang indefinitely.
"""

from __future__ import annotations

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
import urllib.error
import urllib.parse
import urllib.request
import urllib.robotparser
from typing import Any, Iterator

from .errors import RemoteServiceError

__all__ = [
    "HttpClient",
    "HttpResponse",
    "DEFAULT_USER_AGENT",
    "RETRYABLE_STATUS_CODES",
    "cache_dir",
]

log = logging.getLogger(__name__)

RETRYABLE_STATUS_CODES: frozenset[int] = frozenset({408, 425, 429, 500, 502, 503, 504})

DEFAULT_TIMEOUT_SECS: float = 60.0
DEFAULT_MAX_RETRIES: int = 4
DEFAULT_BACKOFF_BASE_SECS: float = 1.0
DEFAULT_BACKOFF_MAX_SECS: float = 30.0
DEFAULT_JITTER_SECS: float = 0.4
# Hard ceiling on the wall-clock time one logical request may consume,
# including every retry and backoff sleep.
DEFAULT_DEADLINE_SECS: float = 120.0
DEFAULT_CHARSET: str = "utf-8"
PROJECT_NAME: str = "agskills"

_FALLBACK_USER_AGENT = (
    "agskills/2.0 (+https://github.com/; academic drug-discovery pipeline)"
)

# A descriptive, contactable User-Agent is required courtesy for public
# scientific APIs and is what their rate-limit policies key on.
DEFAULT_USER_AGENT: str = (
    os.environ.get("AGSKILLS_USER_AGENT")
    or os.environ.get("SCIENCE_SKILLS_USER_AGENT")
    or _FALLBACK_USER_AGENT
)

_THROTTLE_STATUS_RE = re.compile(
    r"(\w[\w ]*?) status:\s*(Green|Yellow|Red|Black)\s*\((\d+)%\)"
)
_THROTTLE_BACKPRESSURE: dict[str, float] = {
    "Green": 0.0, "Yellow": 1.0, "Red": 5.0, "Black": 30.0,
}

# ---------------------------------------------------------------------------
# Cross-platform file locking
# ---------------------------------------------------------------------------

try:  # POSIX
    import fcntl

    def _lock_file(fh) -> None:
        fcntl.flock(fh, fcntl.LOCK_EX)

    def _unlock_file(fh) -> None:
        with contextlib.suppress(OSError):
            fcntl.flock(fh, fcntl.LOCK_UN)

except ImportError:  # Windows
    import msvcrt

    def _lock_file(fh) -> None:
        # Best effort: a contended lock raises rather than blocking forever.
        for _ in range(50):
            try:
                msvcrt.locking(fh.fileno(), msvcrt.LK_NBLCK, 1)
                return
            except OSError:
                time.sleep(0.02)
        # Proceed unlocked rather than deadlock; the limiter is advisory.

    def _unlock_file(fh) -> None:
        with contextlib.suppress(OSError):
            msvcrt.locking(fh.fileno(), msvcrt.LK_UNLCK, 1)


def cache_dir() -> str:
    """Return a writable per-user cache directory, creating it if needed.

    Honours ``AGSKILLS_CACHE_DIR``, then falls back to the platform
    temporary directory. Never returns a hardcoded POSIX path.
    """
    base = os.environ.get("AGSKILLS_CACHE_DIR") or os.path.join(
        tempfile.gettempdir(), PROJECT_NAME
    )
    try:
        os.makedirs(base, exist_ok=True)
    except OSError:
        return tempfile.gettempdir()
    return base


# ---------------------------------------------------------------------------
# Rate limiting
# ---------------------------------------------------------------------------


class _RateLimiter:
    """Enforce a minimum interval between requests to one host.

    State lives in a small file under :func:`cache_dir` so that separate
    skill invocations running concurrently still respect one shared budget.
    The file stores a wall-clock ``time.time()`` value; the lock is released
    before sleeping so concurrent processes interleave rather than queue.
    """

    def __init__(self, hostname: str, qps: float):
        if qps <= 0:
            raise ValueError("qps must be positive")
        self.hostname = hostname
        self.min_interval = 1.0 / qps
        safe_host = re.sub(r"[^A-Za-z0-9._-]", "_", hostname or "unknown")
        self.lock_file = os.path.join(cache_dir(), f"{safe_host}.ratelimit")

    def _claim_slot(self) -> float:
        """Reserve the next slot; return how long the caller should sleep."""
        now = time.time()
        try:
            with open(self.lock_file, "a+", encoding="ascii") as fh:
                _lock_file(fh)
                try:
                    fh.seek(0)
                    raw = fh.read().strip()
                    try:
                        last = float(raw) if raw else 0.0
                    except ValueError:
                        last = 0.0
                    # A timestamp in the future by more than the interval, or
                    # far in the past, means a stale/corrupt file: reset.
                    if last > now + 3600 or last < 0:
                        last = 0.0
                    slot = max(now, last + self.min_interval)
                    fh.seek(0)
                    fh.truncate()
                    fh.write(f"{slot:.6f}")
                    fh.flush()
                finally:
                    _unlock_file(fh)
        except OSError:
            # Cache unavailable: degrade to a simple local delay.
            return self.min_interval
        return max(0.0, slot - now)

    def wait(self, min_sleep: float = 0.0) -> None:
        """Block until the next request to this host is permitted."""
        delay = max(self._claim_slot(), min_sleep)
        if delay > 0:
            time.sleep(delay)


# ---------------------------------------------------------------------------
# robots.txt
# ---------------------------------------------------------------------------


class _RobotsChecker:
    """Lazy robots.txt checker with a working on-disk cache.

    Per RFC 9309 sec. 2.4 an unreachable robots.txt means "allow", so every
    failure path here fails open.
    """

    _CACHE_TTL_SECS: float = 86_400.0

    def __init__(self, base_url: str, user_agent: str, timeout: float = 10.0,
                 cache_directory: str | None = None):
        self._user_agent = user_agent
        self._timeout = timeout
        parsed = urllib.parse.urlparse(base_url)
        self._robots_url = f"{parsed.scheme}://{parsed.netloc}/robots.txt"
        host = re.sub(r"[^A-Za-z0-9._-]", "_", parsed.hostname or "unknown")
        self._cache_file = os.path.join(cache_directory or cache_dir(), f"{host}.robots")
        self._parser: urllib.robotparser.RobotFileParser | None = None
        self._loaded = False

    def _read_cache(self) -> str | None:
        try:
            if time.time() - os.stat(self._cache_file).st_mtime > self._CACHE_TTL_SECS:
                return None
            with open(self._cache_file, "r", encoding="utf-8") as fh:
                return fh.read()
        except OSError:
            return None

    def _write_cache(self, content: str) -> None:
        try:
            with open(self._cache_file, "w", encoding="utf-8") as fh:
                fh.write(content)
        except OSError:
            log.debug("could not cache robots.txt at %s", self._cache_file)

    def _load(self) -> None:
        if self._loaded:
            return
        self._loaded = True
        body = self._read_cache()
        if body is None:
            try:
                req = urllib.request.Request(
                    self._robots_url, headers={"User-Agent": self._user_agent}
                )
                with urllib.request.urlopen(req, timeout=self._timeout) as resp:
                    body = resp.read().decode("utf-8", errors="replace")
            except Exception:
                body = ""  # fail open
            self._write_cache(body)
        parser = urllib.robotparser.RobotFileParser()
        parser.parse(body.splitlines())
        self._parser = parser

    def is_allowed(self, url: str) -> bool:
        self._load()
        if self._parser is None:
            return True
        try:
            return self._parser.can_fetch(self._user_agent, url)
        except Exception:
            return True


# ---------------------------------------------------------------------------
# Responses
# ---------------------------------------------------------------------------


def _maybe_decompress(response):
    encoding = response.headers.get("Content-Encoding", "").lower()
    if encoding in ("gzip", "x-gzip"):
        return gzip.GzipFile(fileobj=response)
    return response


class HttpResponse:
    """A fully-read HTTP response."""

    __slots__ = ("data", "status_code", "headers", "url", "encoding")

    def __init__(self, data: bytes, status_code: int, headers: dict[str, str],
                 url: str, encoding: str | None = None):
        self.data = data
        self.status_code = status_code
        self.headers = headers
        self.url = url
        self.encoding = encoding or DEFAULT_CHARSET

    def json(self) -> Any:
        try:
            return json.loads(self.data.decode(self.encoding, errors="replace"))
        except json.JSONDecodeError as exc:
            preview = self.data[:200].decode("utf-8", errors="replace")
            raise RemoteServiceError(
                f"{self.url} returned {len(self.data)} bytes that are not valid "
                f"JSON: {exc}",
                hint=f"First bytes of the response: {preview!r}",
            ) from exc

    @property
    def text(self) -> str:
        return self.data.decode(self.encoding, errors="replace")

    def __repr__(self) -> str:
        return (f"HttpResponse(status={self.status_code}, url={self.url!r}, "
                f"size={len(self.data)})")


def _parse_retry_after(headers) -> float | None:
    value = headers.get("Retry-After")
    if value is None:
        return None
    try:
        return max(0.0, float(value))
    except (TypeError, ValueError):
        pass
    try:
        when = email.utils.parsedate_to_datetime(value)
        if when.tzinfo is None:
            when = when.replace(tzinfo=datetime.timezone.utc)
        return max(0.0, (when - datetime.datetime.now(datetime.timezone.utc)).total_seconds())
    except (TypeError, ValueError, OverflowError):
        return None


def _parse_throttle_control(headers) -> float:
    value = headers.get("X-Throttling-Control")
    if not value:
        return 0.0
    return max(
        (_THROTTLE_BACKPRESSURE.get(m.group(2), 0.0)
         for m in _THROTTLE_STATUS_RE.finditer(value)),
        default=0.0,
    )


# ---------------------------------------------------------------------------
# Client
# ---------------------------------------------------------------------------


class HttpClient:
    """A rate-limited, retrying HTTP client scoped to one base URL.

    Example:
        >>> client = HttpClient("https://www.ebi.ac.uk/chembl/api/data", qps=5)
        >>> data = client.fetch_json("molecule/CHEMBL25.json")  # doctest: +SKIP
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
        deadline: float = DEFAULT_DEADLINE_SECS,
        user_agent: str = "",
        retryable_status_codes: frozenset[int] = RETRYABLE_STATUS_CODES,
        respect_robots: bool = False,
    ):
        """Create a client scoped to one host.

        ``respect_robots`` defaults to **False**, which is a deliberate and
        narrow choice. robots.txt (RFC 9309) is a directive for automated
        *crawlers* that discover and index content by following links. These
        clients do none of that: they retrieve one named record from a
        documented public REST API because a user asked for it, which is the
        access pattern the providers publish these APIs for. No mainstream
        HTTP client library - ``requests``, ``httpx``, ``urllib3`` - and no
        scientific API client consults robots.txt either.

        The distinction is not academic here. The AlphaFold Database's
        robots.txt contains ``Disallow: /api/``, aimed at crawlers, while
        ``/api/prediction/{accession}`` is the endpoint its own
        documentation tells users to call to locate a model. Honouring the
        crawl directive makes the documented API unusable.

        Set ``respect_robots=True`` for anything that genuinely crawls:
        following links, enumerating pages, or fetching URLs the user did
        not name. Rate limiting, a descriptive and contactable User-Agent,
        and server-directed backoff via ``Retry-After`` apply either way,
        and those are what protect the service.
        """
        parsed = urllib.parse.urlparse(base_url)
        if not parsed.scheme or not parsed.netloc:
            raise ValueError(f"base_url must be absolute, got {base_url!r}")
        self.base_url = base_url
        self.hostname = parsed.hostname or ""
        self.max_retries = max(0, max_retries)
        self.timeout = timeout
        self.backoff_base = backoff_base
        self.backoff_max = backoff_max
        self.jitter = jitter
        self.deadline = deadline
        # Resolved at call time so tests can monkeypatch the env var.
        self.user_agent = user_agent or DEFAULT_USER_AGENT
        self.retryable_status_codes = retryable_status_codes
        self.default_headers = dict(default_headers or {})
        self._limiter = _RateLimiter(self.hostname, qps=qps)
        self._next_min_sleep = 0.0
        self._robots = (
            _RobotsChecker(base_url, self.user_agent) if respect_robots else None
        )

    # -- internals ---------------------------------------------------------

    def _compute_backoff(self, attempt: int, retry_after: float | None = None) -> float:
        delay = self.backoff_base * (2 ** attempt)
        if retry_after is not None:
            delay = max(delay, retry_after)
        delay = min(delay, self.backoff_max)
        if self.jitter > 0:
            delay += random.uniform(0, self.jitter)
        return delay

    def _resolve_url(self, url: str) -> str:
        if "://" not in url:
            base = self.base_url if self.base_url.endswith("/") else self.base_url + "/"
            return urllib.parse.urljoin(base, url.lstrip("/"))
        if urllib.parse.urlparse(url).hostname != self.hostname:
            raise ValueError(
                f"URL {url!r} is not on this client's host {self.hostname!r}"
            )
        return url

    def _build_request(self, url, method, headers, data, json_body):
        merged = {"User-Agent": self.user_agent, "Accept-Encoding": "gzip"}
        merged.update(self.default_headers)
        if headers:
            merged.update(headers)
        body = data
        if json_body is not None:
            body = json.dumps(json_body).encode("utf-8")
            merged.setdefault("Content-Type", "application/json")
        return urllib.request.Request(url, data=body, headers=merged, method=method)

    @contextlib.contextmanager
    def _open_stream(self, url, method, headers, data, json_body, timeout,
                     max_retries=None) -> Iterator[http.client.HTTPResponse]:
        if data is not None and json_body is not None:
            raise ValueError("pass either 'data' or 'json_body', not both")

        url = self._resolve_url(url)
        if self._robots and not self._robots.is_allowed(url):
            raise RemoteServiceError(
                f"robots.txt for {self.hostname} disallows fetching {url}"
            )

        effective_timeout = self.timeout if timeout is None else timeout
        attempts = self.max_retries if max_retries is None else max(0, max_retries)
        started = time.monotonic()
        last_error: str = "no attempt was made"
        next_min_sleep = 0.0

        for attempt in range(attempts + 1):
            if time.monotonic() - started > self.deadline:
                raise RemoteServiceError(
                    f"Gave up on {url} after {self.deadline:.0f}s "
                    f"({attempt} attempt(s)). Last error: {last_error}",
                    hint="The service may be down. Retry later, or use an "
                         "offline subcommand if one is available.",
                )

            self._limiter.wait(min_sleep=max(next_min_sleep, self._next_min_sleep))
            self._next_min_sleep = 0.0
            request = self._build_request(url, method, headers, data, json_body)

            try:
                response = urllib.request.urlopen(request, timeout=effective_timeout)
            except urllib.error.HTTPError as exc:
                status = exc.code
                body = exc.read()
                retry_after = _parse_retry_after(exc.headers)
                exc.close()
                last_error = f"HTTP {status}"
                if status in self.retryable_status_codes and attempt < attempts:
                    next_min_sleep = self._compute_backoff(attempt, retry_after)
                    log.info("%s returned HTTP %d; retrying in >=%.1fs (%d/%d)",
                             url, status, next_min_sleep, attempt + 1, attempts)
                    continue
                hint = None
                if status == 403:
                    hint = ("The server may be rejecting the User-Agent. Set "
                            "AGSKILLS_USER_AGENT to a descriptive string that "
                            "includes a contact address.")
                elif status == 404:
                    hint = "The resource does not exist; check the identifier."
                raise RemoteServiceError(
                    f"HTTP {status} fetching {url}: "
                    f"{body[:300].decode('utf-8', errors='replace')}",
                    hint=hint,
                ) from exc
            except (urllib.error.URLError, OSError, http.client.HTTPException) as exc:
                last_error = f"{type(exc).__name__}: {exc}"
                if attempt < attempts:
                    next_min_sleep = self._compute_backoff(attempt)
                    log.info("network error for %s (%s); retrying in >=%.1fs (%d/%d)",
                             url, exc, next_min_sleep, attempt + 1, attempts)
                    continue
                raise RemoteServiceError(
                    f"Could not reach {url}: {exc}",
                    hint="Check network connectivity and any proxy settings.",
                ) from exc

            throttle = _parse_throttle_control(response.headers)
            if throttle > 0:
                self._next_min_sleep = throttle
            try:
                yield response
            finally:
                response.close()
            return

        raise RemoteServiceError(
            f"Exhausted {attempts} retries for {url}. Last error: {last_error}"
        )

    # -- public API --------------------------------------------------------

    def fetch(self, url: str, *, method: str = "GET",
              headers: dict[str, str] | None = None, data: bytes | None = None,
              json_body: Any | None = None, timeout: float | None = None,
              max_retries: int | None = None) -> HttpResponse:
        """Perform a request and return the fully-read response."""
        with self._open_stream(url, method, headers, data, json_body, timeout,
                               max_retries) as resp:
            body = _maybe_decompress(resp).read()
            charset = resp.headers.get_content_charset() or DEFAULT_CHARSET
            return HttpResponse(body, resp.status, dict(resp.headers),
                                resp.url, charset)

    def fetch_json(self, url: str, **kwargs) -> Any:
        """Fetch *url* and parse the body as JSON."""
        headers = dict(kwargs.pop("headers", None) or {})
        headers.setdefault("Accept", "application/json")
        return self.fetch(url, headers=headers, **kwargs).json()

    def fetch_text(self, url: str, **kwargs) -> str:
        """Fetch *url* and return the decoded body."""
        return self.fetch(url, **kwargs).text

    def fetch_bytes(self, url: str, **kwargs) -> bytes:
        """Fetch *url* and return the raw body."""
        return self.fetch(url, **kwargs).data

    def stream_lines(self, url: str, **kwargs) -> Iterator[str]:
        """Yield the response body line by line without buffering it all."""
        method = kwargs.pop("method", "GET")
        headers = kwargs.pop("headers", None)
        data = kwargs.pop("data", None)
        json_body = kwargs.pop("json_body", None)
        timeout = kwargs.pop("timeout", None)
        max_retries = kwargs.pop("max_retries", None)
        with self._open_stream(url, method, headers, data, json_body, timeout,
                               max_retries) as resp:
            stream = _maybe_decompress(resp)
            charset = resp.headers.get_content_charset() or DEFAULT_CHARSET
            for raw in stream:
                yield raw.decode(charset, errors="replace").rstrip("\r\n")
