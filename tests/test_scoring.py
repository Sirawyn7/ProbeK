import pandas as pd

from probek.models import BlastHit, ClassifiedHit
from probek.scoring import rank_and_select

_DUMMY_HIT = BlastHit(
    sseqid="NC_1", pident=100.0, length=20, qstart=1, qend=20,
    sstart=100, send=119, evalue=1e-10, qlen=20,
)


def _classified(*tiers: str) -> list[ClassifiedHit]:
    return [ClassifiedHit(hit=_DUMMY_HIT, tier=t, locus_label=None) for t in tiers]


def _hit(tier: str, locus_label: str | None = None) -> ClassifiedHit:
    return ClassifiedHit(hit=_DUMMY_HIT, tier=tier, locus_label=locus_label)


def _df(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows)


def test_exon_hit_excludes_regardless_of_quality_and_hervk():
    df = _df(
        [
            {"name": "risky", "sequence": "AAA", "quality": 99.0},
            {"name": "clean", "sequence": "BBB", "quality": 1.0},
        ]
    )
    classified = {
        # High quality and heavy HERV-K coverage, but a real off-target -- must be excluded.
        "AAA": [_hit("A"), _hit("A"), _hit("A"), _hit("A"), _hit("A"), _hit("C", "GENEX")],
        "BBB": [],
    }
    result = rank_and_select("t1", df, classified, top_n=1)
    assert result.selected_names == {"clean"}
    ranks = {r["name"]: r["rank"] for r in result.ranked_rows}
    assert ranks["clean"] == 1
    assert ranks["risky"] == 2


def test_hervk_tiebreak_only_on_exact_quality_tie_fewer_wins():
    df = _df(
        [
            {"name": "few_hervk", "sequence": "AAA", "quality": 50.0},
            {"name": "many_hervk", "sequence": "BBB", "quality": 50.0},
        ]
    )
    classified = {
        "AAA": [_hit("A")],
        "BBB": [_hit("A"), _hit("A"), _hit("A")],
    }
    result = rank_and_select("t1", df, classified, top_n=12)
    order = [r["name"] for r in result.ranked_rows]
    assert order == ["few_hervk", "many_hervk"]


def test_intron_tiebreak_only_when_quality_and_hervk_both_tied():
    df = _df(
        [
            {"name": "few_intron", "sequence": "AAA", "quality": 50.0},
            {"name": "many_intron", "sequence": "BBB", "quality": 50.0},
        ]
    )
    classified = {
        "AAA": [_hit("A"), _hit("A"), _hit("B", "GENEY")],
        "BBB": [_hit("A"), _hit("A"), _hit("B", "GENEY"), _hit("B", "GENEZ"), _hit("B", "GENEZ")],
    }
    result = rank_and_select("t1", df, classified, top_n=12)
    order = [r["name"] for r in result.ranked_rows]
    assert order == ["few_intron", "many_intron"]


def test_outside_gene_hits_never_affect_rank():
    # Identical on every key that matters (quality, HERV-K, intron); only
    # outside_gene_hits differ. If it were used, "fewer_outside" (listed
    # second) would win; since it's ignored, original row order wins instead.
    df = _df(
        [
            {"name": "many_outside", "sequence": "AAA", "quality": 50.0},
            {"name": "fewer_outside", "sequence": "BBB", "quality": 50.0},
        ]
    )
    classified = {
        "AAA": [_hit("B", None), _hit("B", None), _hit("B", None)],
        "BBB": [],
    }
    result = rank_and_select("t1", df, classified, top_n=12)
    order = [r["name"] for r in result.ranked_rows]
    assert order == ["many_outside", "fewer_outside"]


def test_backfill_from_excluded_pool_orders_by_fewest_exon_hits_first():
    df = _df(
        [
            {"name": "e1", "sequence": "E1", "quality": 20.0},
            {"name": "e2", "sequence": "E2", "quality": 10.0},
            {"name": "x1", "sequence": "X1", "quality": 99.0},
            {"name": "x2", "sequence": "X2", "quality": 50.0},
            {"name": "x3", "sequence": "X3", "quality": 40.0},
        ]
    )
    classified = {
        "E1": [],
        "E2": [],
        "X1": [_hit("C", "G1"), _hit("C", "G2")],  # 2 exon hits
        "X2": [_hit("C", "G1")],  # 1 exon hit
        "X3": [_hit("C", "G1")],  # 1 exon hit
    }
    result = rank_and_select("t1", df, classified, top_n=4)

    # 5 total candidates >= top_n=4, so this isn't a "short" target -- it's
    # short on *eligible* candidates specifically, a different condition.
    assert result.short is False

    order = [r["name"] for r in result.ranked_rows]
    # eligible (by quality) first, then excluded (fewest exon hits, then quality)
    assert order == ["e1", "e2", "x2", "x3", "x1"]
    assert result.selected_names == {"e1", "e2", "x2", "x3"}
    ranks = {r["name"]: r["rank"] for r in result.ranked_rows}
    assert ranks == {"e1": 1, "e2": 2, "x2": 3, "x3": 4, "x1": 5}


def test_fully_excluded_target_backfills_entire_selection():
    df = _df([{"name": f"p{i}", "sequence": f"SEQ{i}", "quality": float(i)} for i in range(3)])
    classified = {f"SEQ{i}": [_hit("C", "GENEX")] for i in range(3)}
    result = rank_and_select("t1", df, classified, top_n=12)
    assert result.short is True
    assert result.selected_names == {"p0", "p1", "p2"}


def test_top_n_truncation_and_selection():
    df = _df([{"name": f"p{i}", "sequence": f"SEQ{i}", "quality": float(i)} for i in range(20)])
    classified = {f"SEQ{i}": _classified() for i in range(20)}
    result = rank_and_select("target1", df, classified, top_n=12)
    assert len(result.selected_names) == 12
    assert not result.short
    # highest-quality candidates (19 down to 8) should be selected
    assert result.selected_names == {f"p{i}" for i in range(8, 20)}


def test_short_target_flag_when_fewer_than_top_n():
    df = _df([{"name": f"p{i}", "sequence": f"SEQ{i}", "quality": float(i)} for i in range(5)])
    classified = {f"SEQ{i}": _classified() for i in range(5)}
    result = rank_and_select("target1", df, classified, top_n=12)
    assert result.short is True
    assert len(result.selected_names) == 5  # all candidates selected, none backfilled


def test_rank_assigned_in_order():
    df = _df([{"name": "p1", "sequence": "AAA", "quality": 1.0}, {"name": "p2", "sequence": "BBB", "quality": 2.0}])
    classified = {"AAA": _classified(), "BBB": _classified()}
    result = rank_and_select("target1", df, classified, top_n=12)
    ranks = {r["name"]: r["rank"] for r in result.ranked_rows}
    assert ranks["p2"] == 1  # higher quality ranks first
    assert ranks["p1"] == 2
