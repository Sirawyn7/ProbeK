from probek.chrom_map import ChromMap, parse_assembly_report


def test_parse_assembly_report_drops_na_rows(fixtures_dir):
    df = parse_assembly_report(fixtures_dir / "sample_assembly_report.txt")
    assert set(df["ucsc_name"]) == {"chr1", "chr2"}
    assert "NW_009646194.1" not in set(df["refseq_accn"])


def test_chrom_map_bidirectional(fixtures_dir):
    df = parse_assembly_report(fixtures_dir / "sample_assembly_report.txt")
    cmap = ChromMap.from_dataframe(df)
    assert cmap.to_refseq("chr1") == "NC_000001.11"
    assert cmap.to_ucsc("NC_000002.12") == "chr2"
    assert cmap.to_refseq("chrUnknown") is None
