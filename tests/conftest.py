"""Pytest configuration and shared fixtures for the dnastorage test suite.

The ``test_fi.py`` tests expect ``test_dna.txt`` and ``test_rate.csv`` to exist
in the current working directory.  The ``fi_test_data`` fixture (applied
automatically to every test in this directory) generates those files in a
temporary directory and changes into it so that the tests can locate them with
their original hard-coded relative paths.
"""

import os
import random

import pytest


_STRAND_LEN = 100      # length of each synthetic DNA strand
_NUM_STRANDS = 100     # number of strands written to test_dna.txt
_NUCLEOTIDES = "ACGT"

_OVERALL_ERROR_RATE = 0.01   # 1 % error rate per position
_DEL_RATIO          = 0.333  # fraction of errors that are deletions
_INS_RATIO          = 0.333  # fraction of errors that are insertions
_SUB_RATIO          = 0.334  # fraction of errors that are substitutions


@pytest.fixture(autouse=True)
def fi_test_data(tmp_path, monkeypatch):
    """Create synthetic test-data files and ``chdir`` to *tmp_path*.

    ``test_dna.txt``  – one DNA strand per line (100 strands, each 100 nt).
    ``test_rate.csv`` – per-position error rates used by
                        ``FiTestSuite.test_error_distribution``.

    Using ``monkeypatch.chdir`` ensures that the original working directory is
    restored after each test, even on failure.
    """
    rng = random.Random(42)

    # --- test_dna.txt --------------------------------------------------------
    strands = [
        "".join(rng.choice(_NUCLEOTIDES) for _ in range(_STRAND_LEN))
        for _ in range(_NUM_STRANDS)
    ]
    (tmp_path / "test_dna.txt").write_text("\n".join(strands) + "\n")

    # --- test_rate.csv -------------------------------------------------------
    # Four rows; each row has a label followed by _STRAND_LEN float values.
    # Values at positions 0.._STRAND_LEN-1 are uniform; the inject_distribution
    # method only reads positions p1 .. (strand_len - p2 - 1) from each row.
    def _row(label, value):
        vals = ",".join(str(value) for _ in range(_STRAND_LEN))
        return f"{label},{vals}"

    csv_lines = [
        _row("Overall Error", _OVERALL_ERROR_RATE),
        _row("Del/Error",     _DEL_RATIO),
        _row("Ins/Error",     _INS_RATIO),
        _row("Sub/Error",     _SUB_RATIO),
    ]
    (tmp_path / "test_rate.csv").write_text("\n".join(csv_lines) + "\n")

    # Change into the temp directory so tests can open files by name.
    monkeypatch.chdir(tmp_path)
