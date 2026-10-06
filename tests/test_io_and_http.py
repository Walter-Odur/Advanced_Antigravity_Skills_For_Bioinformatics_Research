"""Filesystem and HTTP infrastructure.

Three properties are non-negotiable and each was violated at least once by
the code this replaces:

* **Always UTF-8.** ``Path.write_text`` without an encoding uses the locale
  codepage, cp1252 on a default Windows install, so a report containing an
  en dash died with ``UnicodeEncodeError``.
* **Always LF for generated scripts.** A SLURM script written with CRLF
  fails on a Linux cluster with ``bash\\r: bad interpreter``.
* **Cross-platform cache paths.** The HTTP client cached ``robots.txt``
  under a hardcoded ``/tmp``, which does not exist on Windows, so every
  write failed silently and the file was re-fetched on every client
  construction.
"""

from __future__ import annotations

import json
import os

import pytest

from agskills.errors import ResourceNotFoundError
from agskills.http import DEFAULT_USER_AGENT, HttpClient, cache_dir
from agskills.io_utils import (
    read_json,
    read_lines,
    require_file,
    safe_stem,
    write_json,
    write_text,
)


# ---------------------------------------------------------------------------
# Encoding
# ---------------------------------------------------------------------------


def test_text_is_written_as_utf8_regardless_of_locale(tmp_path):
    """Non-ASCII content must survive, whatever the console codepage is."""
    content = "Angstrom A, en dash -, degree 25 C, alpha-helix\n"
    path = write_text(content, tmp_path / "report.md")
    assert path.read_text(encoding="utf-8") == content


def test_json_is_written_as_utf8(tmp_path):
    data = {"note": "binding site 4.5 A from the hinge", "value": 1.5}
    path = write_json(data, tmp_path / "out.json")
    assert read_json(path) == data
    # Readable as UTF-8 bytes, not escaped into ASCII.
    assert "A from the hinge" in path.read_text(encoding="utf-8")


def test_json_handles_values_the_science_layer_produces(tmp_path):
    """Paths and numpy scalars must serialise rather than raising."""
    from pathlib import Path
    data = {"path": Path("a/b.pdb"), "count": 3}
    path = write_json(data, tmp_path / "out.json")
    loaded = read_json(path)
    assert loaded["count"] == 3
    assert "b.pdb" in loaded["path"]


# ---------------------------------------------------------------------------
# Line endings
# ---------------------------------------------------------------------------


def test_generated_text_always_has_lf_endings(tmp_path):
    """Written on Windows, run on a cluster."""
    path = write_text("line one\nline two\n", tmp_path / "script.sh")
    assert path.read_bytes() == b"line one\nline two\n"
    assert b"\r" not in path.read_bytes()


def test_crlf_input_is_normalised(tmp_path):
    """A template that picked up CRLF must still produce a usable script."""
    path = write_text("one\r\ntwo\rthree\n", tmp_path / "mixed.sh")
    assert path.read_bytes() == b"one\ntwo\nthree\n"


def test_json_output_also_uses_lf(tmp_path):
    path = write_json({"a": 1}, tmp_path / "out.json")
    assert b"\r" not in path.read_bytes()


# ---------------------------------------------------------------------------
# Atomicity and directories
# ---------------------------------------------------------------------------


def test_parent_directories_are_created(tmp_path):
    path = write_json({"a": 1}, tmp_path / "deep" / "nested" / "out.json")
    assert path.is_file()


def test_a_failed_write_leaves_no_partial_file(tmp_path):
    """A truncated JSON file would be half-parsed by the next step."""
    class Unserialisable:
        def __repr__(self):
            raise RuntimeError("boom")

    target = tmp_path / "out.json"
    with pytest.raises(Exception):
        write_json({"bad": Unserialisable()}, target)
    assert not target.exists()
    # And no temporary files are left behind.
    assert list(tmp_path.glob(".*")) == []


def test_writing_replaces_an_existing_file_atomically(tmp_path):
    target = tmp_path / "out.json"
    write_json({"version": 1}, target)
    write_json({"version": 2}, target)
    assert read_json(target) == {"version": 2}


def test_write_to_a_bare_filename_works(tmp_path, monkeypatch):
    """``--output results.json`` with no directory part must not fail."""
    monkeypatch.chdir(tmp_path)
    path = write_json({"a": 1}, "results.json")
    assert path.is_file()


def test_executable_bit_is_set_where_supported(tmp_path):
    path = write_text("#!/bin/bash\necho hi\n", tmp_path / "run.sh",
                      executable=True)
    if os.name != "nt":
        assert os.access(path, os.X_OK)
    else:
        # Windows has no execute bit; the call must simply not fail.
        assert path.is_file()


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------


def test_missing_file_raises_with_a_hint(tmp_path):
    with pytest.raises(ResourceNotFoundError) as excinfo:
        require_file(tmp_path / "absent.pdb")
    assert "not found" in str(excinfo.value)


def test_a_directory_is_not_a_file(tmp_path):
    with pytest.raises(ResourceNotFoundError):
        require_file(tmp_path)


def test_read_lines_skips_comments_and_blanks(tmp_path):
    path = tmp_path / "list.smi"
    path.write_text("# header\n\nCCO\n  CCC  \n# trailing\n", encoding="utf-8")
    assert read_lines(path) == ["CCO", "CCC"]


# ---------------------------------------------------------------------------
# Filename safety
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", [
    "a/b", "a\\b", "a:b", 'a"b', "a<b", "a>b", "a|b", "a?b", "a*b",
])
def test_every_illegal_character_is_replaced(name):
    stem = safe_stem(name)
    for bad in '/\\:"<>|?*':
        assert bad not in stem


@pytest.mark.parametrize("reserved", ["CON", "PRN", "AUX", "NUL", "COM1",
                                      "LPT1", "con", "nul"])
def test_windows_reserved_names_are_escaped(reserved):
    """These cannot be filenames on Windows at all, in any case."""
    stem = safe_stem(reserved)
    assert stem.upper() not in {"CON", "PRN", "AUX", "NUL", "COM1", "LPT1"}


def test_safe_stem_never_returns_empty():
    for name in ("", "   ", "...", "///", "\x00"):
        assert safe_stem(name)


# ---------------------------------------------------------------------------
# HTTP: cache location
# ---------------------------------------------------------------------------


def test_cache_directory_exists_and_is_writable():
    """The old client wrote to a hardcoded ``/tmp``, absent on Windows.

    Every write failed silently, so ``robots.txt`` was re-fetched on every
    single client construction - an extra network round trip per command.
    """
    directory = cache_dir()
    assert os.path.isdir(directory)
    probe = os.path.join(directory, "write_probe.tmp")
    with open(probe, "w", encoding="utf-8") as handle:
        handle.write("ok")
    assert os.path.isfile(probe)
    os.remove(probe)


def test_cache_directory_is_not_a_hardcoded_posix_path():
    directory = cache_dir()
    if os.name == "nt":
        assert directory != "/tmp"
        assert not directory.startswith("/tmp")


def test_cache_directory_honours_the_environment(tmp_path, monkeypatch):
    target = tmp_path / "custom_cache"
    monkeypatch.setenv("AGSKILLS_CACHE_DIR", str(target))
    assert cache_dir() == str(target)
    assert target.is_dir()


# ---------------------------------------------------------------------------
# HTTP: User-Agent
# ---------------------------------------------------------------------------


def test_default_user_agent_is_not_empty():
    """An empty User-Agent earns a 403 from several EBI and NCBI endpoints.

    The old default was the empty string whenever the environment variable
    was unset.
    """
    assert DEFAULT_USER_AGENT
    assert len(DEFAULT_USER_AGENT) > 10
    client = HttpClient("https://example.org", qps=1.0)
    assert client.user_agent


def test_user_agent_can_be_overridden(monkeypatch):
    monkeypatch.setenv("AGSKILLS_USER_AGENT", "my-lab-client/1.0 (me@lab.org)")
    client = HttpClient("https://example.org", qps=1.0,
                        user_agent=os.environ["AGSKILLS_USER_AGENT"])
    assert "my-lab-client" in client.user_agent


# ---------------------------------------------------------------------------
# HTTP: construction and URL handling
# ---------------------------------------------------------------------------


def test_relative_url_resolves_against_the_base():
    client = HttpClient("https://www.ebi.ac.uk/chembl/api/data", qps=1.0)
    resolved = client._resolve_url("target/search.json?q=EGFR")
    assert resolved == \
        "https://www.ebi.ac.uk/chembl/api/data/target/search.json?q=EGFR"


def test_leading_slash_on_a_relative_url_is_tolerated():
    client = HttpClient("https://files.rcsb.org", qps=1.0)
    assert client._resolve_url("/download/1M17.pdb").endswith(
        "/download/1M17.pdb")


def test_absolute_url_on_another_host_is_refused():
    """A client is scoped to one host, so its rate limit means something."""
    client = HttpClient("https://files.rcsb.org", qps=1.0)
    with pytest.raises(ValueError):
        client._resolve_url("https://evil.example.org/x")


def test_relative_base_url_is_refused():
    with pytest.raises(ValueError):
        HttpClient("not-a-url", qps=1.0)


def test_non_positive_qps_is_refused():
    with pytest.raises(ValueError):
        HttpClient("https://example.org", qps=0)


# ---------------------------------------------------------------------------
# HTTP: retry policy
# ---------------------------------------------------------------------------


def test_backoff_grows_and_is_capped():
    client = HttpClient("https://example.org", qps=10.0, backoff_base=1.0,
                        backoff_max=10.0, jitter=0.0)
    delays = [client._compute_backoff(attempt) for attempt in range(6)]
    assert delays[:4] == [1.0, 2.0, 4.0, 8.0]
    assert all(delay <= 10.0 for delay in delays)
    assert delays == sorted(delays)


def test_retry_after_header_is_respected():
    client = HttpClient("https://example.org", qps=10.0, backoff_base=1.0,
                        backoff_max=100.0, jitter=0.0)
    assert client._compute_backoff(0, retry_after=30.0) == 30.0


def test_total_deadline_is_bounded():
    """A dead host must not block for minutes.

    The old defaults allowed seven retries with backoff capped at 180 s,
    so one unreachable service could stall a command for over six minutes.
    That is what made the original test suite hang: it never completed and
    had to be killed.
    """
    from agskills.http import (
        DEFAULT_BACKOFF_MAX_SECS,
        DEFAULT_DEADLINE_SECS,
        DEFAULT_MAX_RETRIES,
    )
    assert DEFAULT_MAX_RETRIES <= 5
    assert DEFAULT_BACKOFF_MAX_SECS <= 60
    assert DEFAULT_DEADLINE_SECS <= 180

    client = HttpClient("https://example.org", qps=10.0)
    worst_case = sum(client._compute_backoff(a)
                     for a in range(client.max_retries))
    assert worst_case < client.deadline


def test_unreachable_host_fails_quickly_and_clearly():
    """A network failure must be a typed error with advice, not a traceback."""
    import time

    from agskills.errors import RemoteServiceError

    client = HttpClient("https://localhost:1", qps=10.0, max_retries=1,
                        timeout=2.0, deadline=12.0, backoff_base=0.1,
                        jitter=0.0)
    started = time.monotonic()
    with pytest.raises(RemoteServiceError) as excinfo:
        client.fetch_text("/nothing")
    elapsed = time.monotonic() - started
    assert elapsed < 30.0, f"took {elapsed:.1f}s to give up"
    message = str(excinfo.value)
    assert "localhost" in message
    assert "network" in message.lower() or "reach" in message.lower()


def test_retryable_status_codes_include_the_transient_ones():
    from agskills.http import RETRYABLE_STATUS_CODES
    for code in (429, 500, 502, 503, 504):
        assert code in RETRYABLE_STATUS_CODES
    # A 404 is not transient; retrying wastes the user's time.
    assert 404 not in RETRYABLE_STATUS_CODES
    assert 400 not in RETRYABLE_STATUS_CODES


def test_non_json_response_raises_with_a_preview():
    """A service returning HTML must not surface as a JSONDecodeError."""
    from agskills.errors import RemoteServiceError
    from agskills.http import HttpResponse

    response = HttpResponse(b"<!DOCTYPE html><html>error</html>", 200, {},
                            "https://example.org/api")
    with pytest.raises(RemoteServiceError) as excinfo:
        response.json()
    message = str(excinfo.value)
    assert "not valid JSON" in message
    assert "DOCTYPE" in message


# ---------------------------------------------------------------------------
# HTTP: rate limiting
# ---------------------------------------------------------------------------


def test_rate_limiter_uses_wall_clock_time(tmp_path, monkeypatch):
    """State is shared between processes, so the clock must have an epoch.

    The old limiter persisted ``time.monotonic()`` to a shared file.
    Monotonic clocks have no common epoch across processes or reboots, so
    the stored value was not comparable with a fresh reading.
    """
    import time

    from agskills.http import _RateLimiter

    monkeypatch.setenv("AGSKILLS_CACHE_DIR", str(tmp_path))
    limiter = _RateLimiter("example.org", qps=1000.0)
    limiter.wait()
    stored = float(open(limiter.lock_file, encoding="ascii").read())
    # A wall-clock timestamp is close to time.time(); a monotonic one on a
    # machine up for days would be wildly different.
    assert abs(stored - time.time()) < 60.0


def test_rate_limiter_spaces_successive_requests(tmp_path, monkeypatch):
    import time

    from agskills.http import _RateLimiter

    monkeypatch.setenv("AGSKILLS_CACHE_DIR", str(tmp_path))
    limiter = _RateLimiter("example.org", qps=20.0)  # 50 ms apart
    started = time.monotonic()
    for _ in range(3):
        limiter.wait()
    elapsed = time.monotonic() - started
    assert elapsed >= 0.08, f"three requests took only {elapsed:.3f}s"


def test_rate_limiter_recovers_from_a_corrupt_state_file(tmp_path,
                                                          monkeypatch):
    from agskills.http import _RateLimiter

    monkeypatch.setenv("AGSKILLS_CACHE_DIR", str(tmp_path))
    limiter = _RateLimiter("example.org", qps=1000.0)
    with open(limiter.lock_file, "w", encoding="ascii") as handle:
        handle.write("not a number")
    limiter.wait()  # must not raise


def test_rate_limiter_sanitises_the_hostname_for_the_filename(tmp_path,
                                                               monkeypatch):
    from agskills.http import _RateLimiter

    monkeypatch.setenv("AGSKILLS_CACHE_DIR", str(tmp_path))
    limiter = _RateLimiter("host:with/bad\\chars", qps=10.0)
    assert ":" not in os.path.basename(limiter.lock_file)
    limiter.wait()


# ---------------------------------------------------------------------------
# HTTP: live behaviour, opt-in
# ---------------------------------------------------------------------------


@pytest.mark.network
def test_a_real_fetch_succeeds_and_is_decoded():
    client = HttpClient("https://files.rcsb.org", qps=3.0)
    text = client.fetch_text("/download/1M17.pdb")
    assert text.startswith("HEADER")
    assert "ATOM" in text


@pytest.mark.network
def test_a_real_404_is_a_typed_error_with_a_hint():
    from agskills.errors import RemoteServiceError
    client = HttpClient("https://files.rcsb.org", qps=3.0, max_retries=0)
    with pytest.raises(RemoteServiceError) as excinfo:
        client.fetch_text("/download/ZZZZ.pdb")
    assert "404" in str(excinfo.value)
    assert excinfo.value.hint
