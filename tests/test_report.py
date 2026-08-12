import pandas as pd

from probek.models import BlastHit, ClassifiedHit, TargetResult
from probek.report import (
    _reorder_columns,
    enrich_target_result,
    feature_counts,
    flagged_genes,
    format_off_target_loci,
    off_target_risk_label,
)

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
    assert row["off_target_risk"] == "Low"
    assert row["flagged_genes"] == ""
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


def test_flagged_genes_deduplicates_and_ignores_non_exon_tiers():
    classified = [
        ClassifiedHit(hit=_DUMMY_HIT, tier="C", locus_label="GM2A"),
        ClassifiedHit(hit=_DUMMY_HIT, tier="C", locus_label="GM2A"),  # duplicate, same gene
        ClassifiedHit(hit=_DUMMY_HIT, tier="C", locus_label="TP53"),
        ClassifiedHit(hit=_DUMMY_HIT, tier="A", locus_label="HERVK-int"),  # not an exon hit
        ClassifiedHit(hit=_DUMMY_HIT, tier="B", locus_label="ARMC3"),  # intron, not exon
    ]
    assert flagged_genes(classified) == "GM2A, TP53"


def test_flagged_genes_empty_when_no_exon_hits():
    classified = [ClassifiedHit(hit=_DUMMY_HIT, tier="A", locus_label="HERVK-int")]
    assert flagged_genes(classified) == ""


def test_off_target_risk_label_thresholds():
    assert off_target_risk_label(0) == "Low"
    assert off_target_risk_label(1) == "Moderate"
    assert off_target_risk_label(2) == "Moderate"
    assert off_target_risk_label(3) == "High"
    assert off_target_risk_label(10) == "High"


def test_reorder_columns_puts_priority_columns_first_and_loci_last():
    df = pd.DataFrame(
        [
            {
                "length": 21,
                "sequence": "ACGT",
                "off_target_loci": "some detail",
                "name": "p1",
                "quality": 50.0,
                "exon_hits": 0,
                "selected": True,
                "rank": 1,
                "off_target_risk": "Low",
            }
        ]
    )
    reordered = _reorder_columns(df)
    columns = list(reordered.columns)
    assert columns[-1] == "off_target_loci"
    assert columns.index("name") < columns.index("length")
    assert columns.index("selected") < columns.index("sequence")
    assert columns.index("off_target_risk") < columns.index("quality")


def test_reorder_columns_handles_missing_priority_columns_gracefully():
    # final_selection.csv-only columns ("target") absent from a plain audit df
    df = pd.DataFrame([{"name": "p1", "off_target_loci": "x", "length": 21}])
    reordered = _reorder_columns(df)
    assert list(reordered.columns) == ["name", "length", "off_target_loci"]
