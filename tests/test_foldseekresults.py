import pytest
import os

# import shutil
from plaseek.tools.foldseek import (
    retrieve_foldseek_results,
    write_merged_m8file,
)


@pytest.fixture(scope="function")
def fixture1():
    result_file, database = retrieve_foldseek_results(
        "tests/inputfiles/AF-P07676-F1-model_v4.pdb", mode="3diaa"
    )
    yield result_file
    os.remove(result_file)
    os.remove("alis_merged.m8")


def test_write_merged_m8file(fixture1):
    write_merged_m8file(fixture1)
    assert os.path.exists("alis_merged.m8")
    with open("alis_merged.m8") as f:
        assert "AF-A0A4D4JGC3-F1-model_v4" in f.read()
        f.seek(0)
        assert "MGYP003606827913.pdb.gz" in f.read()
