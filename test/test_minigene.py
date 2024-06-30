from pathlib import Path
from process_test_data import process_test_data
import pandas as pd
import pytest
from difflib import SequenceMatcher

# Run these tests from the project root directory with `python -m pytest`
# Or, if you're using VS Code, use the Testing tab

@pytest.fixture(scope='session')
def results():
    HG38FOLDER = Path("resources")
    INPUTFOLDER = Path("test")
    return process_test_data(HG38FOLDER, INPUTFOLDER)

def test_data_exists(results):
    assert isinstance(results['ref'], pd.DataFrame)
    assert isinstance(results['var'], pd.DataFrame)
    assert isinstance(results['both'], pd.DataFrame)

def test_n_records(results):
    assert len(results['ref']) == len(results['var'])
    assert len(results['both'] == len(results['ref']) + len(results['var']))
    assert len(results['ref'].columns) + 1 == len(results['both'].columns)

def test_all_minigenes_93(results):
    for refvar in ['ref', 'var']:
        counted = results[refvar].copy()
        counted['minigene_length'] = counted['minigene'].str.len() == 93
        assert counted['minigene_length'].all()

def test_inserted_diff(results):
    """
    For every insertion result, compare the
    ref and var minigenes and confirm that
     - there is a differing substring within the var sequence
     - that substring is equal to the variant codon
     - the ref and var minigenes are otherwise identical
    """
    df = results['both']
    df = df.loc[df['variant_class'] == 'insertion']
    for index, row in df.iterrows():
        ref_gene = row['minigene_ref']
        var_gene = row['minigene_var']
        matcher = SequenceMatcher(a=ref_gene, b=var_gene, autojunk=False)
        blocks = matcher.get_matching_blocks()
        # Should be 2-3 blocks (2 if the diff is the first character, three otherwise)
        if len(blocks) == 1:
            raise Exception(f"Ref and Var minigenes at [{index}] are identical")
        # Check with Kelvin - Looking at the output CSV, why doesn't the var codon match what ends up in the minigene?
        # print(f"Expected {row['variant_class']} ({row['codon_ref']} - {row['codon_var']}")
        # print("ref: " + ref_gene)
        # print("var: " + var_gene)
        print(f"{index} has {len(blocks)} blocks")