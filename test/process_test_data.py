from pathlib import Path
from pyfaidx import Fasta
import pandas as pd

from scripts.generate_sequences import generate_sequences
from scripts._utils import extract_exon_info

def process_test_data(hg38_folder: Path, input_folder: Path) -> dict[str, pd.DataFrame]:
    """
    Returns an object containing
        "ref": dataframe containing the reference minigenes
        "var": dataframe containing the variant minigenes
        "both": ref and var dataframes merged on "minigene_id"
    
    The "ref" and "var" dataframes are straight out of the library (they end up as the CSV files).
    The "both" frame is created here for convenience.
    """
    # import reference and data
    print("Importing reference data...")
    fasta_file = hg38_folder / "gdc" / "GRCh38.d1.vd1" / "fasta" / "GRCh38.d1.vd1.fa"
    gtf_file = hg38_folder / "annotation" / "gencode.v33.primary_assembly.annotation.gtf"
    refgen = Fasta(filename=fasta_file)
    exon_info = extract_exon_info(gtf_file)
    input_sample = [("test", input_folder / f"test_hc_vep.txt")]

    # Assemble the sequences for the test sample
    print("Generating sequences...")
    ref_df, var_df = generate_sequences(refgen, exon_info, input_sample)

    merged = without_refvar(ref_df).merge(
        without_refvar(var_df), 
        on="minigene_id", 
        suffixes=("", "_drop"))
    merged = merged.rename({ 
        "minigene": "minigene_ref", 
        "minigene_drop": "minigene_var" }, axis=1)
    dupes = [col for col in merged.columns if col.endswith("_drop")]
    merged = merged.drop(dupes, axis=1)

    return {
        "ref": ref_df, "var": var_df,
        "both": merged
    }

# Removes the _ref and _var from the end of the minigene_id
def without_refvar(df: pd.DataFrame):
    result = df.copy(True)
    result['minigene_id'] = result['minigene_id'].str[:-len('_ref')]
    return result