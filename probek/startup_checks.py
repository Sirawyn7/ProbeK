from __future__ import annotations

import logging
from pathlib import Path

from .exceptions import InputDataError

logger = logging.getLogger(__name__)

DEFAULT_INPUT_DIR = "input_csvs"


def check_input_dir(input_path: Path, explicit: bool) -> bool:
    """Validates the input path. Returns True if there are CSVs to process.

    If this is the default `input_csvs/` dir and it's empty, prints a
    reminder and returns False so the caller can exit(0) cleanly rather than
    erroring — there's simply nothing to do yet. An explicitly-passed path
    with no CSVs, or a path that doesn't exist at all, is a genuine error.
    """
    if not input_path.exists():
        raise InputDataError(f"Input path does not exist: {input_path}")

    if input_path.is_dir():
        csvs = list(input_path.glob("*.csv"))
        if not csvs:
            if not explicit and input_path.name == DEFAULT_INPUT_DIR:
                logger.info(
                    "No CSVs found in %s/ yet.\n"
                    "Drop your eFISHent-format CSVs there (one per gene target, or "
                    "combined) and re-run. See %s/README.md for the expected format.",
                    input_path,
                    input_path,
                )
                return False
            raise InputDataError(f"No CSV files found in: {input_path}")

    return True


__all__ = ["check_input_dir", "DEFAULT_INPUT_DIR"]
