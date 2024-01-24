import pytest
import tempfile

# import shutil
from plaseek.plaseek import parse_tsvfile, filtering_by_pident


@pytest.fixture
def open_foldseekhits():
    foldseekhits = parse_tsvfile("tests/inputfiles/AF-P07676-F1-model_v4.tsv")
    yield foldseekhits


def test_parse_tsvdile(open_foldseekhits):
    """test for parse_tsvfile() function."""
    assert len(open_foldseekhits) == 650
    assert open_foldseekhits[0].id == "AF-P07676-F1-model_v4"
    assert len(open_foldseekhits[0].seq) == 382


def test_filtering_by_pident():
    """test for filtering_by_pident() function."""
    tsvfile = "tests/inputfiles/test_intermediate.tsv"
    output = tempfile.NamedTemporaryFile(suffix=".txt", delete=False).name
    filtering_by_pident(
        tsvfile,
        pident_threshold=98.0,
    )
    assert open(output).read() == open("tests/inputfiles/test_output.txt").read()
