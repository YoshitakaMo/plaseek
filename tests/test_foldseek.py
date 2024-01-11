# %%
import pytest
from repseek import parse_tsvfile


@pytest.fixture
def open_foldseekhits():
    foldseekhits = parse_tsvfile("tests/inputfiles/AF-P07676-F1-model_v4.tsv")
    yield foldseekhits


# %%
# foldseekhits = parse_tsvfile("inputfiles/AF-P07676-F1-model_v4.tsv")

# %%
# len(foldseekhits[0].seq)


# %%
def test_parse_tsvdile(open_foldseekhits):
    """test for parse_tsvfile() function."""
    assert len(open_foldseekhits) == 650
    assert open_foldseekhits[0].id == "AF-P07676-F1-model_v4"
    assert len(open_foldseekhits[0].seq) == 382
