import pytest
import tempfile
from repseek import parse_tsvfile, collect_pident_plasmid


@pytest.fixture
def open_foldseekhits():
    foldseekhits = parse_tsvfile("tests/inputfiles/AF-P07676-F1-model_v4.tsv")
    yield foldseekhits


def test_parse_tsvdile(open_foldseekhits):
    """test for parse_tsvfile() function."""
    assert len(open_foldseekhits) == 650
    assert open_foldseekhits[0].id == "AF-P07676-F1-model_v4"
    assert len(open_foldseekhits[0].seq) == 382
    # TODO: check whether the file contains
    # "qaccver","saccver", "pident", "length", "evalue", "bitscore" colomns.
    assert open_foldseekhits[0].saccver == "NC_000913.3"


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
    """test for run_foldseek() function."""
    fastafile = "tests/inputfiles/rk2_TrfA2.fasta"
    dbtype = "pdb"
    run_foldseek()
