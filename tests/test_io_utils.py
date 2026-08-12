import pandas as pd

from probek.io_utils import dedup_sequences, derive_target_label, drop_failed, load_input_csvs


def test_derive_target_label_from_common_prefix():
    names = ["pNRV101_gag-1", "pNRV101_gag-2", "pNRV101_gag-16"]
    assert derive_target_label(names, "some_filename") == "pNRV101_gag"


def test_derive_target_label_falls_back_to_filename_on_disagreement():
    names = ["pNRV101_gag-1", "pNRV102_pro-1"]
    assert derive_target_label(names, "mixed_file") == "mixed_file"


def test_derive_target_label_falls_back_when_no_trailing_index():
    names = ["probe_alpha", "probe_beta"]
    assert derive_target_label(names, "fallback_name") == "fallback_name"


def test_load_input_csvs_derives_target_and_drops_nothing_yet(fixtures_dir):
    frames = load_input_csvs(fixtures_dir / "sample_efishent.csv")
    assert set(frames.keys()) == {"pNRV101_gag"}
    assert len(frames["pNRV101_gag"]) == 3  # FAIL row still present pre-filter


def test_drop_failed_removes_fail_rows_only():
    df = pd.DataFrame(
        {
            "name": ["a", "b", "c"],
            "recommendation": ["FAIL", "", "FLAG(low_complexity)"],
        }
    )
    result = drop_failed(df)
    assert list(result["name"]) == ["b", "c"]


def test_dedup_sequences_merges_identical_sequences_across_targets():
    df_gag = pd.DataFrame(
        {"name": ["pNRV101_gag-1"], "sequence": ["ACGTACGTACGTACGTACGT"]}
    )
    df_pro = pd.DataFrame(
        {"name": ["pNRV102_pro-16"], "sequence": ["ACGTACGTACGTACGTACGT"]}
    )
    groups = dedup_sequences({"pNRV101_gag": df_gag, "pNRV102_pro": df_pro})
    assert len(groups) == 1
    [group] = groups.values()
    assert len(group.occurrences) == 2
    targets = {occ.target for occ in group.occurrences}
    assert targets == {"pNRV101_gag", "pNRV102_pro"}


def test_dedup_sequences_is_case_insensitive_key():
    df = pd.DataFrame({"name": ["a", "b"], "sequence": ["acgt", "ACGT"]})
    groups = dedup_sequences({"target": df})
    assert len(groups) == 1
