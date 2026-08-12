from __future__ import annotations

import hashlib
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

from .models import BlastHit

_SCHEMA = """
CREATE TABLE IF NOT EXISTS blast_cache (
    sequence              TEXT PRIMARY KEY,
    seq_length            INTEGER NOT NULL,
    hits_json             TEXT NOT NULL,
    reference_db_version  TEXT NOT NULL,
    blast_params_hash     TEXT NOT NULL,
    blastn_version        TEXT NOT NULL,
    created_at             TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_blast_cache_refver ON blast_cache(reference_db_version);
"""


def open_cache(db_path: Path) -> sqlite3.Connection:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)
    conn.executescript(_SCHEMA)
    conn.commit()
    return conn


def blast_params_hash(params: dict) -> str:
    encoded = json.dumps(params, sort_keys=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def get_cached_hits(
    conn: sqlite3.Connection,
    sequence: str,
    reference_db_version: str,
    params_hash: str,
) -> list[BlastHit] | None:
    """Return cached hits, or None on a cache miss (absent, or stale relative
    to the current reference data version / BLAST params)."""
    row = conn.execute(
        "SELECT hits_json, reference_db_version, blast_params_hash "
        "FROM blast_cache WHERE sequence = ?",
        (sequence,),
    ).fetchone()
    if row is None:
        return None
    hits_json, cached_ref_version, cached_params_hash = row
    if cached_ref_version != reference_db_version or cached_params_hash != params_hash:
        return None
    return [BlastHit.from_dict(d) for d in json.loads(hits_json)]


def store_hits(
    conn: sqlite3.Connection,
    sequence: str,
    seq_length: int,
    hits: list[BlastHit],
    reference_db_version: str,
    params_hash: str,
    blastn_version: str,
) -> None:
    conn.execute(
        "INSERT OR REPLACE INTO blast_cache "
        "(sequence, seq_length, hits_json, reference_db_version, blast_params_hash, "
        " blastn_version, created_at) VALUES (?, ?, ?, ?, ?, ?, ?)",
        (
            sequence,
            seq_length,
            json.dumps([h.to_dict() for h in hits]),
            reference_db_version,
            params_hash,
            blastn_version,
            datetime.now(timezone.utc).isoformat(),
        ),
    )
    conn.commit()


def partition_by_cache(
    sequences: list[str],
    conn: sqlite3.Connection,
    reference_db_version: str,
    params_hash: str,
    force: bool = False,
) -> tuple[list[str], dict[str, list[BlastHit]]]:
    """Split sequences into (needs_blast, cached_hits_by_sequence)."""
    if force:
        return list(sequences), {}
    to_blast: list[str] = []
    cached: dict[str, list[BlastHit]] = {}
    for seq in sequences:
        hits = get_cached_hits(conn, seq, reference_db_version, params_hash)
        if hits is None:
            to_blast.append(seq)
        else:
            cached[seq] = hits
    return to_blast, cached
