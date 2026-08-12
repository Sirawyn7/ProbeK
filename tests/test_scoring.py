import pandas as pd

from probek.models import BlastHit, ClassifiedHit
from probek.scoring import rank_and_select

_DUMMY_HIT = BlastHit(
    sseqid="NC_1", pident=100.0, length=20, qstart=1, qend=20,
    sstart=100, send=119, evalue=1e-10, qlen=20,
)


def _classified(*tiers: str) -> list[ClassifiedHit]:
    return [ClassifiedHit(hit=_DUMMY_HIT, tier=t, locus_label=None) for t in tiers]


def _df(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows)


def test_sort_order_tier_c_then_tier_a_then_tier_b_then_quality():
    df = _df(
        [
            {"name": "low_c_high_a", "sequence": "AAA", "quality": 5.0},
            {"name": "high_c", "sequence": "BBB", "quality": 9.0},
            {"name": "tiebreak_lower_quality", "sequence": "CCC", "quality": 1.0},
            {"name": "tiebreak_higher_quality", "sequence": "DDD", "quality": 8.0},
        ]
    )
    classified = {
        "AAA": _classified("A", "A", "A"),  # tier_c=0, tier_a=3
        "BBB": _classified("C", "C"),  # tier_c=2 -> worst
        "CCC": _classified(),  # tier_c=0, tier_a=0 -> ties with DDD except quality
        "DDD": _classified(),  # tier_c=0, tier_a=0, higher quality -> should rank above CCC
    }
    result = rank_and_select("target1", df, classified, top_n=12)
    order = [r["name"] for r in result.ranked_rows]
    # low_c_high_a: tier_c=0 (best), beats DDD/CCC which also have tier_c=0 but tier_a=0
    assert order[0] == "low_c_high_a"
    # DDD (quality 8) ranks above CCC (quality 1) as the tiebreaker
    assert order.index("tiebreak_higher_quality") < order.index("tiebreak_lower_quality")
    # high_c (tier_c=2) is worst
    assert order[-1] == "high_c"


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
