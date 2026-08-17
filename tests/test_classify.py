import pytest

from probek.classify import classify_hits, feature_counts
from probek.interval_backend import IntervalFeature, build_index
from probek.models import BlastHit


@pytest.fixture(autouse=True)
def force_custom_backend(monkeypatch):
    # Tier-precedence logic is backend-agnostic; pin to the always-available
    # pure pandas/numpy fallback so these tests don't depend on whether
    # pyranges/pyranges1 happen to be installed.
    monkeypatch.setenv("PROBEK_INTERVAL_BACKEND", "custom")


def _hit(sseqid: str, start: int, end: int) -> BlastHit:
    return BlastHit(
        sseqid=sseqid,
        pident=100.0,
        length=end - start + 1,
        qstart=1,
        qend=end - start + 1,
        sstart=start,
        send=end,
        evalue=1e-10,
        qlen=end - start + 1,
    )


@pytest.fixture
def indexes():
    ervk_index = build_index([IntervalFeature(chrom="NC_1", start=1000, end=2000, name="LTR5_Hs")])
    exon_index = build_index(
        [
            IntervalFeature(chrom="NC_1", start=1500, end=1800, name="GENE1"),
            IntervalFeature(chrom="NC_1", start=5000, end=5100, name="GENE2"),
        ]
    )
    gene_index = build_index(
        [
            IntervalFeature(chrom="NC_1", start=4900, end=5200, name="GENE2"),
            IntervalFeature(chrom="NC_1", start=8000, end=8500, name="GENE3"),
        ]
    )
    return ervk_index, exon_index, gene_index


def test_ervk_overlap_wins_over_exon_overlap(indexes):
    # Hit overlaps both the ERVK interval (1000-2000) and GENE1's exon (1500-1800)
    hit = _hit("NC_1", 1600, 1700)
    [classified] = classify_hits([hit], *indexes)
    assert classified.tier == "A"
    assert classified.locus_label == "LTR5_Hs"


def test_exon_overlap_without_ervk_is_tier_c(indexes):
    hit = _hit("NC_1", 5050, 5060)
    [classified] = classify_hits([hit], *indexes)
    assert classified.tier == "C"
    assert classified.locus_label == "GENE2"


def test_intron_only_overlap_is_tier_b_with_nearest_gene(indexes):
    hit = _hit("NC_1", 8100, 8200)
    [classified] = classify_hits([hit], *indexes)
    assert classified.tier == "B"
    assert classified.locus_label == "GENE3"


def test_no_overlap_at_all_is_tier_b_with_no_gene(indexes):
    hit = _hit("NC_1", 20000, 20010)
    [classified] = classify_hits([hit], *indexes)
    assert classified.tier == "B"
    assert classified.locus_label is None


def test_boundary_touching_counts_as_overlap(indexes):
    # ERVK interval ends at 2000; a hit starting exactly there should still overlap.
    hit = _hit("NC_1", 2000, 2001)
    [classified] = classify_hits([hit], *indexes)
    assert classified.tier == "A"


def test_reverse_strand_hit_classified_via_normalized_span(indexes):
    reverse_hit = BlastHit(
        sseqid="NC_1", pident=100.0, length=101, qstart=1, qend=101,
        sstart=1700, send=1600, evalue=1e-10, qlen=101,
    )
    [classified] = classify_hits([reverse_hit], *indexes)
    assert classified.tier == "A"


def test_feature_counts_tallies_all_four_buckets(indexes):
    hits = [
        _hit("NC_1", 1600, 1700),  # HERV-K family
        _hit("NC_1", 5050, 5060),  # exon
        _hit("NC_1", 8100, 8200),  # intron
        _hit("NC_1", 20000, 20010),  # outside any gene
    ]
    classified = classify_hits(hits, *indexes)
    assert feature_counts(classified) == {"hervk": 1, "exon": 1, "intron": 1, "outside_gene": 1}


def test_feature_counts_empty():
    assert feature_counts([]) == {"hervk": 0, "exon": 0, "intron": 0, "outside_gene": 0}


def test_classify_hits_empty_list(indexes):
    assert classify_hits([], *indexes) == []
