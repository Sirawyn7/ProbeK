"""Shared streamed-download helpers with a live progress bar (tqdm)."""

from __future__ import annotations

import hashlib
from pathlib import Path

import requests
from tqdm import tqdm

from .exceptions import DownloadError

_CHUNK_SIZE = 1 << 20  # 1 MiB


def _write_response_to_file(resp: requests.Response, dest: Path, bar: tqdm) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "wb") as f:
        for chunk in resp.iter_content(chunk_size=_CHUNK_SIZE):
            if not chunk:
                continue
            f.write(chunk)
            bar.update(len(chunk))


def download_with_progress(url: str, dest: Path, desc: str | None = None) -> None:
    """Streams `url` to `dest`, showing a live progress bar against the
    response's Content-Length. Falls back to an indeterminate bar (bytes
    transferred and elapsed time, no percentage/ETA) if the server doesn't
    report a length."""
    with requests.get(url, stream=True, timeout=60) as resp:
        resp.raise_for_status()
        total = int(resp.headers.get("content-length") or 0) or None
        with tqdm(
            total=total,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            desc=desc or dest.name,
        ) as bar:
            _write_response_to_file(resp, dest, bar)


def download_volumes_with_progress(
    volumes: list[tuple[str, Path]], total_bytes: int, desc: str
) -> None:
    """Streams multiple URLs to their destinations under a single combined
    progress bar, against a pre-known total size — for cases like a
    multi-volume BLAST database where the exact total is published in a
    metadata file rather than discoverable per-file up front."""
    with tqdm(total=total_bytes, unit="B", unit_scale=True, unit_divisor=1024, desc=desc) as bar:
        for url, dest in volumes:
            with requests.get(url, stream=True, timeout=60) as resp:
                resp.raise_for_status()
                _write_response_to_file(resp, dest, bar)


def verify_md5(path: Path, md5_url: str) -> None:
    """Verifies `path` against the MD5 checksum published at `md5_url`
    (NCBI's standard "<hash>  <filename>" convention)."""
    expected = requests.get(md5_url, timeout=30).text.strip().split()[0]
    actual = hashlib.md5(path.read_bytes()).hexdigest()
    if actual != expected:
        raise DownloadError(
            f"{path.name} failed checksum verification "
            f"(expected {expected}, got {actual}) — try again."
        )
