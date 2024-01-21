import pytest
import tempfile

# import shutil
# import os
# from pathlib import Path
from repseek.repseek import parse_tsvfile, collect_pident_plasmid  # type: ignore


@pytest.fixture
def open_foldseekhits():
    foldseekhits = parse_tsvfile("tests/inputfiles/AF-P07676-F1-model_v4.tsv")
    yield foldseekhits


def test_parse_tsvdile(open_foldseekhits):
    """test for parse_tsvfile() function."""
    assert len(open_foldseekhits) == 650
    assert open_foldseekhits[0].id == "AF-P07676-F1-model_v4"
    assert len(open_foldseekhits[0].seq) == 382


def test_collect_pident_plasmid():
    """test for collect_pident_plasmid() function."""
    tsvfile = "tests/inputfiles/test_intermediate.tsv"
    output = tempfile.NamedTemporaryFile(suffix=".txt", delete=False).name
    collect_pident_plasmid(
        tsvfile,
        output,
        pident_threshold=98.0,
    )
    assert open(output).read() == open("tests/inputfiles/test_output.txt").read()


def test_run_foldseek():
    pass


#     """test for run_foldseek() function.
#     This test requires foldseek binary and foldseek database.
#     yayoi11-14 nodes have foldseek binary and database.
#     """
#     input = Path("tests/inputfiles/AF-P07676-F1-model_v4.pdb")
#     foldseek_binary_path = shutil.which("foldseek")
#     foldseek_binary_path = "/home/apps/foldseek/20231027/bin/foldseek"
#     foldseek_db_path = os.getenv("FOLDSEEKDB")
#     foldseek_db_path = "/scr/foldseek"
#     foldseek_tsvfile = f"{input.stem}.tsv"
#     run_foldseek(
#         pdbfile=input,
#         foldseek_binary_path=foldseek_binary_path,
#         foldseek_db_path=foldseek_db_path,
#         outtsvfile=foldseek_tsvfile,
#     )
