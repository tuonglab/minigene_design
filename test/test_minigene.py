from pathlib import Path
from process_test_data import process_test_data
import pandas as pd
import pytest

# Run these tests from the project root directory with `python -m pytest`
# Or, if you're using VS Code, use the Testing tab


@pytest.fixture(scope="session")
def results():
    HG38FOLDER = Path(__file__).parent.parent / "resources"
    INPUTFOLDER = Path(__file__).parent.parent / "test"
    return process_test_data(HG38FOLDER, INPUTFOLDER)


def test_data_exists(results):
    assert isinstance(results["ref"], pd.DataFrame)
    assert isinstance(results["var"], pd.DataFrame)
    assert isinstance(results["both"], pd.DataFrame)


def test_ref_var_length(results):
    assert len(results["ref"]) == len(results["var"])
    assert len(results["both"] == len(results["ref"]) + len(results["var"]))
    assert len(results["ref"].columns) + 1 == len(results["both"].columns)
