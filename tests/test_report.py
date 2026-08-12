from probek.models import BlastHit, ClassifiedHit, TargetResult
from probek.report import enrich_target_result, feature_counts, format_off_target_loci

_DUMMY_HIT = BlastHit(
    sseqid="NC_000005.10", pident=100.0, length=21, qstart=1, qend=21,
    sstart=151268868, send=151268889, evalue=1e-10, qlen=21,
)


def _hit_at(sseqid: str, start: int, end: int) -> BlastHit:
    return BlastHit(
        sseqid=sseqid, pident=100.0, length=end - start + 1, qstart=1,
        qend=end - start + 1, sstart=start, send=end, evalue=1e-10, qlen=end - start + 1,
    )


def test_plain_labels_cover_all_four_feature_categories():
    classified = [
        ClassifiedHit(hit=_hit_at("NC_1", 100, 120), tier="A", locus_label="HERVK-int"),
        ClassifiedHit(hit=_hit_at("NC_1", 200, 220), tier="C", locus_label="GM2A"),
        ClassifiedHit(hit=_hit_at("NC_1", 300, 320), tier="B", locus_label="ARMC3"),
        ClassifiedHit(hit=_hit_at("NC_1", 400, 420), tier="B", locus_label=None),
    ]
    text = format_off_target_loci(classified)
    parts = text.split("; ")
    assert parts[0] == "NC_1:100-120 (HERV-K family: HERVK-int)"
    assert parts[1] == "NC_1:200-220 (exon of GM2A)"
    assert parts[2] == "NC_1:300-320 (intron of ARMC3)"
    assert parts[3] == "NC_1:400-420 (outside any gene)"


def test_format_off_target_loci_empty():
    assert format_off_target_loci([]) == ""


def test_feature_counts_tallies_all_four_buckets():
    classified = [
        ClassifiedHit(hit=_DUMMY_HIT, tier="A", locus_label="HERVK-int"),
        ClassifiedHit(hit=_DUMMY_HIT, tier="A", locus_label="LTR5_Hs"),
        ClassifiedHit(hit=_DUMMY_HIT, tier="C", locus_label="GM2A"),
        ClassifiedHit(hit=_DUMMY_HIT, tier="B", locus_label="ARMC3"),
        ClassifiedHit(hit=_DUMMY_HIT, tier="B", locus_label=None),
        ClassifiedHit(hit=_DUMMY_HIT, tier="B", locus_label=None),
    ]
    assert feature_counts(classified) == {"hervk": 2, "exon": 1, "intron": 1, "outside_gene": 2}


def test_feature_counts_empty():
    assert feature_counts([]) == {"hervk": 0, "exon": 0, "intron": 0, "outside_gene": 0}


def test_enrich_target_result_replaces_tier_counts_with_plain_columns():
    result = TargetResult(
        target="t1",
        ranked_rows=[
            {
                "name": "p1",
                "sequence": "acgt",
                "tier_a_count": 3,
                "tier_b_count": 1,
                "tier_c_count": 0,
            }
        ],
        selected_names={"p1"},
        short=False,
    )
    classified_by_sequence = {
        "ACGT": [
            ClassifiedHit(hit=_DUMMY_HIT, tier="A", locus_label="HERVK-int"),
            ClassifiedHit(hit=_DUMMY_HIT, tier="A", locus_label="HERVK-int"),
            ClassifiedHit(hit=_DUMMY_HIT, tier="A", locus_label="HERVK-int"),
            ClassifiedHit(hit=_DUMMY_HIT, tier="B", locus_label=None),
        ]
    }
    enrich_target_result(result, classified_by_sequence)

    row = result.ranked_rows[0]
    assert "tier_a_count" not in row
    assert "tier_b_count" not in row
    assert "tier_c_count" not in row
    assert row["hervk_family_hits"] == 3
    assert row["exon_hits"] == 0
    assert row["intron_hits"] == 0
    assert row["outside_gene_hits"] == 1
    assert row["selected"] is True
    assert "(HERV-K family: HERVK-int)" in row["off_target_loci"]


def test_enrich_target_result_sequence_lookup_is_case_insensitive():
    result = TargetResult(
        target="t1",
        ranked_rows=[{"name": "p1", "sequence": "acgt"}],
        selected_names=set(),
        short=False,
    )
    classified_by_sequence = {"ACGT": [ClassifiedHit(hit=_DUMMY_HIT, tier="C", locus_label="GM2A")]}
    enrich_target_result(result, classified_by_sequence)
    assert result.ranked_rows[0]["exon_hits"] == 1
    assert result.ranked_rows[0]["selected"] is False
