import pytest
import tempfile

# import shutil
from plaseek.plaseek import filtering_m8file, filtering_by_pident


@pytest.fixture
def open_foldseekhits():
    foldseekhits = filtering_m8file("tests/inputfiles/AF-P07676-F1-model_v4.tsv")
    yield foldseekhits


def test_parse_tsvdile(open_foldseekhits):
    """test for filtering_m8file() function."""
    assert len(open_foldseekhits) == 650
    assert open_foldseekhits[0].id == "AF-P07676-F1-model_v4"
    assert len(open_foldseekhits[0].seq) == 382


def test_filtering_by_pident():
    """test for filtering_by_pident() function.
    The 50 sequences listed in the validation file must be included in the output file.
    """
    # validation = "tests/validation/IncP-2_lists.txt"
    tsvfile = "tests/inputfiles/test_intermediate.tsv"
    output = tempfile.NamedTemporaryFile(suffix=".txt", delete=False).name
    filtering_by_pident(
        tsvfile,
        pident_threshold=98.0,
    )
    assert open(output).read() == open("tests/inputfiles/test_output.txt").read()
