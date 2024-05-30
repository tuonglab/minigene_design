from pathlib import Path
from pyfaidx import Fasta
import pandas as pd
from scripts._utils import (
    find,
    read_and_filter,
    create_result_list,
    extract_result,
    get_sequences_indel,
    get_sequences_substitution,
)

def generate_sequences(
        refgen: Fasta, exon_info: dict,
        input_samples: list[tuple[str, Path | str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    
    results = {}
    for sample, file_path in input_samples:
        mut_dict = create_result_list(read_and_filter(file_path))
        for mutation_class in get_sequences:
            mutations = find(mut_dict["variant_class"], mutation_class)
            results[sample] = {}
            if len(mutations) > 0:
                for mut in mutations:
                    results[sample][mut] = {}
                    results[sample][mut]["ref"], results[sample][mut]["var"] = (
                        get_sequences[mutation_class](
                            mut_info=extract_result(mut_dict, mut),
                            exon_info=exon_info,
                            fasta=refgen,
                        )
                    )
                    results[sample][mut]["mut_info"] = extract_result(mut_dict, mut)
    
    return create_minigene_dfs(results)

# Maps mutation class to the function that handles it
get_sequences = {
    "insertion": get_sequences_indel,
    "deletion": get_sequences_indel,
    "SNV": get_sequences_substitution,
    "substitution": get_sequences_substitution
}

def create_minigene_dfs(results_dict):
    records = { 
        "ref": [], "var": []
    }
    for refvar in ["ref", "var"]:
        for sample in results_dict:
            for mut in results_dict[sample]:
                if results_dict[sample][mut][refvar] is not None:
                    for protein in results_dict[sample][mut][refvar]:
                        for frame in results_dict[sample][mut][refvar][protein]:
                            for i, seq in enumerate(
                                results_dict[sample][mut][refvar][protein][frame]
                            ):
                                mut_info = {}
                                for key, record in results_dict[sample][mut][
                                    "mut_info"
                                ].items():
                                    mut_info[key] = record
                                mut_info.update(
                                    {
                                        "minigene": seq,
                                        "minigene_id": str(sample)
                                        + "_"
                                        + str(mut)
                                        + "_"
                                        + protein
                                        + "_"
                                        + str(frame)
                                        + "_"
                                        + str(i)
                                        + "_"
                                        + refvar
                                        ,
                                    }
                                )
                                records[refvar].append(mut_info)

    return (pd.DataFrame(records['ref']), pd.DataFrame(records['var']))