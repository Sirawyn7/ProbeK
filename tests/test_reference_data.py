from probek.reference_data import parse_blastdb_metadata


def test_parse_blastdb_metadata_converts_ftp_to_https_and_sorts_by_filename():
    # Real shape of NCBI's <db>-nucl-metadata.json; note "files" isn't
    # guaranteed sorted (NCBI returns .01 before .00 for human_genome).
    metadata = {
        "files": [
            "ftp://ftp.ncbi.nlm.nih.gov/blast/db/human_genome.01.tar.gz",
            "ftp://ftp.ncbi.nlm.nih.gov/blast/db/human_genome.00.tar.gz",
        ],
        "bytes-total-compressed": 1024270191,
    }
    volumes, total = parse_blastdb_metadata(metadata)
    assert total == 1024270191
    assert volumes == [
        ("https://ftp.ncbi.nlm.nih.gov/blast/db/human_genome.00.tar.gz", "human_genome.00.tar.gz"),
        ("https://ftp.ncbi.nlm.nih.gov/blast/db/human_genome.01.tar.gz", "human_genome.01.tar.gz"),
    ]


def test_parse_blastdb_metadata_single_volume():
    metadata = {
        "files": ["ftp://ftp.ncbi.nlm.nih.gov/blast/db/some_small_db.tar.gz"],
        "bytes-total-compressed": 12345,
    }
    volumes, total = parse_blastdb_metadata(metadata)
    assert total == 12345
    assert volumes == [
        ("https://ftp.ncbi.nlm.nih.gov/blast/db/some_small_db.tar.gz", "some_small_db.tar.gz"),
    ]
