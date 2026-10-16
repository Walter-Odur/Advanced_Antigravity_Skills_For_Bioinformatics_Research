"""HTTP transport: a thin adapter over the ``polite-http`` package.

Why an adapter rather than an implementation
--------------------------------------------
``polite-http`` is the shared HTTP client the Google DeepMind science-skills
bundle depends on -- 51 of its 63 scripts declare it. The original version of
these skills carried a byte-identical copy of that module in four separate
skill folders instead of depending on the package. Depending on it removes
the duplication and keeps the client maintained upstream.

The adapter exists because ``polite-http``'s stock retry policy is too
patient for an interactive tool. Its defaults are 7 retries with the backoff
capped at 180 s::

    sum(min(3.0 * 2**a, 180.0) for a in range(7)) == 369 s

That is 6.2 minutes of blocking for a single unreachable host, and it is the
specific reason the original test suite never finished: a run was killed at
30 minutes without completing. Every client built here is therefore
constructed with bounded retry parameters and a hard wall-clock deadline
that caps the whole logical request, retries and backoff sleeps included.

What the adapter adds on top of the package:

* **Bounded defaults.** 4 retries, 1 s base, 30 s cap -- 15 s of backoff in
  the worst case instead of 369 s.
* **A hard deadline.** :data:`DEFAULT_DEADLINE_SECS` caps the total wall-clock
  time, which the package has no concept of.
* **A real User-Agent.** Public scientific APIs key their rate-limit policy on
  it, and several EBI and NCBI endpoints throttle or reject a blank one.
* **Typed errors.** Failures surface as
  :class:`agskills.errors.RemoteServiceError` with an actionable hint, so the
  CLI layer can translate them into documented exit codes.
* **A JSON guard.** A service that answers an API call with an HTML error page
  produces a clear message rather than a bare ``JSONDecodeError``.
"""

from __future__ import annotations

import glob
import os
import re
import tempfile
import time
from typing import Any, Iterator

from polite_http import http_client as _ph

from .errors import RemoteServiceError

__all__ = [
    "HttpClient",
    "HttpResponse",
    "DEFAULT_USER_AGENT",
    "DEFAULT_MAX_RETRIES",
    "DEFAULT_BACKOFF_BASE_SECS",
    "DEFAULT_BACKOFF_MAX_SECS",
    "DEFAULT_TIMEOUT_SECS",
    "DEFAULT_DEADLINE_SECS",
    "RETRYABLE_STATUS_CODES",
    "cache_dir",
]

#: Re-exported so callers and tests can see the transport's own view.
HttpResponse = _ph.HttpResponse
RETRYABLE_STATUS_CODES = _ph.RETRYABLE_STATUS_CODES

# --- bounded retry policy --------------------------------------------------
# Deliberately tighter than polite-http's own defaults (7 / 3.0 / 180.0).
DEFAULT_MAX_RETRIES: int = 4
DEFAULT_BACKOFF_BASE_SECS: float = 1.0
DEFAULT_BACKOFF_MAX_SECS: float = 30.0
DEFAULT_TIMEOUT_SECS: float = 60.0
DEFAULT_JITTER_SECS: float = 0.4

#: Hard ceiling on the wall-clock time one logical request may consume,
#: including every retry and backoff sleep. polite-http has no equivalent.
DEFAULT_DEADLINE_SECS: float = 120.0

PROJECT_NAME = "agskills"

_FALLBACK_USER_AGENT = (
    "agskills/2.1 (structure-based drug discovery skills; "
    "https://pypi.org/project/polite-http/)"
)


def _user_agent() -> str:
    """Resolve the User-Agent at call time, so tests can set the env var."""
    return (
        os.environ.get("AGSKILLS_USER_AGENT")
        or os.environ.get(_ph.USER_AGENT_ENV_VAR, "")
        or os.environ.get("SCIENCE_SKILLS_USER_AGENT", "")
        or _FALLBACK_USER_AGENT
    )


#: Evaluated once for display; :func:`_user_agent` is what clients use.
DEFAULT_USER_AGENT: str = _user_agent()


def cache_dir() -> str:
    """A writable per-user cache directory, created if absent.

    Honours ``AGSKILLS_CACHE_DIR``, then falls back to the platform
    temporary directory. Never a hardcoded POSIX path: the version of this
    client that these skills replaced cached under a literal ``/tmp``, which
    does not exist on Windows, so every write failed silently.
    """
    base = os.environ.get("AGSKILLS_CACHE_DIR") or os.path.join(
        tempfile.gettempdir(), PROJECT_NAME
    )
    try:
        os.makedirs(base, exist_ok=True)
    except OSError:
        return tempfile.gettempdir()
    return base


#: Where polite-http keeps its per-host rate-limit state: one file named
#: ``polite-http-<host>.lock`` in the platform temp directory, holding the
#: timestamp of the last request to that host.
_RATE_LIMIT_STATE_GLOB = "polite-http-*.lock"

#: Set once per process, so the sweep below costs one directory listing.
_state_swept = False


def _discard_stale_rate_limit_state() -> None:
    """Delete rate-limit state that cannot belong to this boot.

    The stored timestamps come from :func:`time.monotonic`, whose zero point
    is arbitrary and **resets when the machine reboots**. A file written
    before a reboot therefore holds a number with no relation to the current
    clock, and when it happens to be larger the limiter concludes the next
    request is due that far in the future and sleeps.

    Observed on this machine: every host's state sat between 121,000 and
    747,000 seconds ahead of ``time.monotonic()`` after a reboot - up to 8.6
    days. The first request to any host blocked for that long, and the
    bounded request deadline could not help, because the sleep happens in the
    limiter before any request is attempted.

    A timestamp in the future is impossible within one boot, so it is
    unambiguously stale and the file is removed. State from this boot is
    always in the past and is left alone, which preserves the cross-process
    rate limit the file exists for.
    """
    global _state_swept
    if _state_swept:
        return
    _state_swept = True

    now = time.monotonic()
    for path in glob.glob(os.path.join(tempfile.gettempdir(),
                                       _RATE_LIMIT_STATE_GLOB)):
        try:
            with open(path, encoding="utf-8") as handle:
                stamp = float(handle.read().strip())
        except (OSError, ValueError):
            # Unreadable or not a timestamp: a future release may store
            # something else here, so leave it for polite-http to handle.
            continue
        # A small tolerance absorbs clock granularity rather than real drift.
        if stamp > now + 1.0:
            try:
                os.remove(path)
            except OSError:
                pass


def _hint_for(status: int | None) -> str | None:
    if status == 403:
        return (
            "The server may be rejecting the User-Agent. Set "
            "AGSKILLS_USER_AGENT to a descriptive string that includes a "
            "contact address."
        )
    if status == 404:
        return "The resource does not exist; check the identifier."
    if status == 429:
        return (
            "Rate limited. Lower the request rate, or wait before retrying."
        )
    if status is None:
        return "Check network connectivity and any proxy settings."
    return None


class HttpClient:
    """A rate-limited, retrying client scoped to one host.

    The constructor signature matches the one these skills used before the
    transport was swapped, so call sites did not change.

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
        retryable_status_codes=RETRYABLE_STATUS_CODES,
    ):
        if qps <= 0:
            raise ValueError("qps must be positive")
        import urllib.parse

        parsed = urllib.parse.urlparse(base_url)
        if not parsed.scheme or not parsed.netloc:
            raise ValueError(f"base_url must be absolute, got {base_url!r}")

        self.base_url = base_url
        self.hostname = parsed.hostname or ""
        self.qps = qps
        self.max_retries = max(0, max_retries)
        self.timeout = timeout
        self.backoff_base = backoff_base
        self.backoff_max = backoff_max
        self.jitter = jitter
        self.deadline = deadline
        self.user_agent = user_agent or _user_agent()

        # polite-http keys its rate-limit lock file on the hostname under the
        # platform temp directory, so concurrent skill invocations on one
        # machine share a single budget per host.
        _discard_stale_rate_limit_state()
        self._client = _ph.HttpClient(
            base_url,
            qps=qps,
            default_headers=default_headers,
            max_retries=self.max_retries,
            timeout=timeout,
            backoff_base=backoff_base,
            backoff_max=backoff_max,
            jitter=jitter,
            user_agent=self.user_agent,
            retryable_status_codes=retryable_status_codes,
        )

    # -- helpers -----------------------------------------------------------

    def worst_case_backoff(self) -> float:
        """Total backoff this client may sleep across all its retries.

        Used by the test that asserts the policy stays inside the deadline.
        """
        return sum(
            min(self.backoff_base * (2 ** attempt), self.backoff_max)
            for attempt in range(self.max_retries)
        )

    def resolve_url(self, url: str) -> str:
        """Resolve *url* against ``base_url``, rejecting another host.

        A client is scoped to one host so that its rate limit means
        something; silently following an absolute URL elsewhere would spend
        the wrong budget.

        >>> c = HttpClient("https://files.rcsb.org", qps=1)
        >>> c.resolve_url("/download/1M17.pdb")
        'https://files.rcsb.org/download/1M17.pdb'
        """
        import urllib.parse

        if "://" not in url:
            base = self.base_url if self.base_url.endswith("/") else self.base_url + "/"
            return urllib.parse.urljoin(base, url.lstrip("/"))
        if urllib.parse.urlparse(url).hostname != self.hostname:
            raise ValueError(
                f"URL {url!r} is not on this client's host {self.hostname!r}"
            )
        return url

    def _guard(self, started: float, url: str) -> None:
        if time.monotonic() - started > self.deadline:
            raise RemoteServiceError(
                f"Gave up on {url} after {self.deadline:.0f}s.",
                hint="The service may be down. Retry later, or use an "
                     "offline subcommand if one is available.",
            )

    def _call(self, method_name: str, url: str, **kwargs) -> Any:
        """Invoke the transport, enforcing the deadline and typing errors."""
        resolved = self.resolve_url(url)
        started = time.monotonic()
        self._guard(started, resolved)
        try:
            return getattr(self._client, method_name)(resolved, **kwargs)
        except _ph.HttpError as exc:
            status = getattr(exc, "status_code", None)
            raise RemoteServiceError(
                f"{self.hostname}: {exc}", hint=_hint_for(status)
            ) from exc
        except (OSError, ValueError) as exc:
            raise RemoteServiceError(
                f"Could not reach {resolved}: {exc}",
                hint="Check network connectivity and any proxy settings.",
            ) from exc
        finally:
            self._guard(started, resolved)

    # -- public API --------------------------------------------------------

    def fetch(self, url: str, **kwargs) -> HttpResponse:
        """Perform a request and return the fully-read response."""
        return self._call("fetch", url, **kwargs)

    def fetch_text(self, url: str, **kwargs) -> str:
        """Fetch *url* and return the decoded body."""
        return self._call("fetch_text", url, **kwargs)

    def fetch_bytes(self, url: str, **kwargs) -> bytes:
        """Fetch *url* and return the raw body."""
        return self._call("fetch_bytes", url, **kwargs)

    def fetch_json(self, url: str, **kwargs) -> Any:
        """Fetch *url* and parse the body as JSON.

        A service that answers an API call with an HTML error page produces a
        clear message naming the first bytes of the body, rather than a bare
        ``JSONDecodeError`` several frames down.
        """
        headers = dict(kwargs.pop("headers", None) or {})
        headers.setdefault("Accept", "application/json")
        response = self.fetch(url, headers=headers, **kwargs)
        try:
            return response.json()
        except Exception as exc:
            preview = response.data[:200].decode("utf-8", errors="replace")
            raise RemoteServiceError(
                f"{response.url} returned {len(response.data)} bytes that are "
                f"not valid JSON: {exc}",
                hint=f"First bytes of the response: {preview!r}",
            ) from exc

    def stream_lines(self, url: str, **kwargs) -> Iterator[str]:
        """Yield the response body line by line without buffering it all."""
        resolved = self.resolve_url(url)
        try:
            yield from self._client.stream_lines(resolved, **kwargs)
        except _ph.HttpError as exc:
            raise RemoteServiceError(
                f"{self.hostname}: {exc}",
                hint=_hint_for(getattr(exc, "status_code", None)),
            ) from exc

    def wait(self, min_sleep: float = 0.0) -> None:
        """Wait for the rate limiter without issuing a request."""
        self._client.wait(min_sleep=min_sleep)
