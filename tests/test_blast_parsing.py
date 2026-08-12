from probek.blast import filter_hits, parse_outfmt6


def test_parse_outfmt6_groups_by_query(fixtures_dir):
    hits = parse_outfmt6(fixtures_dir / "sample_blast_outfmt6.tsv")
    assert set(hits.keys()) == {"SEQAAA", "SEQBBB"}
    assert len(hits["SEQAAA"]) == 3
    assert len(hits["SEQBBB"]) == 1


def test_coverage_computation(fixtures_dir):
    hits = parse_outfmt6(fixtures_dir / "sample_blast_outfmt6.tsv")
    hit = hits["SEQAAA"][0]
    assert hit.length == 20
    assert hit.qlen == 20
    assert hit.coverage == 1.0


def test_span_normalizes_reverse_strand(fixtures_dir):
    hits = parse_outfmt6(fixtures_dir / "sample_blast_outfmt6.tsv")
    reverse_hit = hits["SEQAAA"][2]
    assert reverse_hit.sstart == 2000
    assert reverse_hit.send == 1983
    assert reverse_hit.span == (1983, 2000)


def test_filter_hits_drops_below_threshold(fixtures_dir):
    hits = parse_outfmt6(fixtures_dir / "sample_blast_outfmt6.tsv")
    genuine = filter_hits(hits["SEQAAA"], min_identity=90.0, min_coverage=90.0)
    # row 2 (85% identity) must be dropped; rows 1 and 3 (100% identity, >=90% coverage) survive
    assert len(genuine) == 2
    assert all(h.pident >= 90.0 for h in genuine)
    assert all(h.coverage * 100 >= 90.0 for h in genuine)


def test_filter_hits_empty_input():
    assert filter_hits([], min_identity=90.0, min_coverage=90.0) == []
