from probek.cache import (
    blast_params_hash,
    get_cached_hits,
    open_cache,
    partition_by_cache,
    store_hits,
)
from probek.models import BlastHit

_HIT = BlastHit(
    sseqid="NC_1", pident=100.0, length=20, qstart=1, qend=20,
    sstart=100, send=119, evalue=1e-10, qlen=20,
)


def _cache(tmp_path):
    return open_cache(tmp_path / "cache.sqlite3")


def test_cache_roundtrip(tmp_path):
    conn = _cache(tmp_path)
    params_hash = blast_params_hash({"task": "blastn-short"})
    store_hits(conn, "ACGT", 4, [_HIT], "refv1", params_hash, "blastn 2.16.0")

    hits = get_cached_hits(conn, "ACGT", "refv1", params_hash)
    assert hits == [_HIT]


def test_cache_miss_when_absent(tmp_path):
    conn = _cache(tmp_path)
    assert get_cached_hits(conn, "NOTCACHED", "refv1", "hash1") is None


def test_cache_miss_on_reference_version_mismatch(tmp_path):
    conn = _cache(tmp_path)
    params_hash = blast_params_hash({"task": "blastn-short"})
    store_hits(conn, "ACGT", 4, [_HIT], "refv1", params_hash, "blastn 2.16.0")

    assert get_cached_hits(conn, "ACGT", "refv2", params_hash) is None


def test_cache_miss_on_params_hash_mismatch(tmp_path):
    conn = _cache(tmp_path)
    params_hash = blast_params_hash({"task": "blastn-short"})
    store_hits(conn, "ACGT", 4, [_HIT], "refv1", params_hash, "blastn 2.16.0")

    other_hash = blast_params_hash({"task": "blastn-short", "word_size": 11})
    assert get_cached_hits(conn, "ACGT", "refv1", other_hash) is None


def test_partition_by_cache_splits_hit_and_miss(tmp_path):
    conn = _cache(tmp_path)
    params_hash = blast_params_hash({"task": "blastn-short"})
    store_hits(conn, "CACHED", 6, [_HIT], "refv1", params_hash, "blastn 2.16.0")

    to_blast, cached = partition_by_cache(
        ["CACHED", "NEW"], conn, "refv1", params_hash, force=False
    )
    assert to_blast == ["NEW"]
    assert cached == {"CACHED": [_HIT]}


def test_partition_by_cache_force_bypasses_everything(tmp_path):
    conn = _cache(tmp_path)
    params_hash = blast_params_hash({"task": "blastn-short"})
    store_hits(conn, "CACHED", 6, [_HIT], "refv1", params_hash, "blastn 2.16.0")

    to_blast, cached = partition_by_cache(
        ["CACHED", "NEW"], conn, "refv1", params_hash, force=True
    )
    assert set(to_blast) == {"CACHED", "NEW"}
    assert cached == {}
