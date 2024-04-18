# From neoepitopes files (class ii) eg. /scratch/project/tcr_neoantigen/results/results_healthyXtumour/neoantigens/sample1/Class_I/sample1_MHC_Class_I_all_epitopes_ccf.tsv
# https://stackoverflow.com/questions/64243197/how-to-query-a-region-in-a-fasta-file-using-biopython
# Get 45bp upstream and downstream of the seq (from neoepitopes) (unique sequence filter)

import pandas as pd
from pyfaidx import Fasta
from collections import defaultdict

# import reference and data
reference = Fasta(
    filename="/scratch/project/tcr_neoantigen/resources/references/hg38/gdc/GRCh38.d1.vd1/fasta/GRCh38.d1.vd1.fa"
)
df = pd.read_csv(
    "/scratch/project/tcr_neoantigen/results/blood_results_171123/neoantigens/bloodXTumour/Class_I/bloodXTumour_MHC_Class_I_all_epitopes_ccf_ref_match.tsv",
    sep="\t",
)
df = df.drop_duplicates(["Stop", "Reference", "Variant"])
# create lists of mutation location and varients
mut_list = list(df["Stop"])
chrom_list = (list(df["Chromosome"]),)
ref_list = list(df["Reference"])
var_list = list(df["Variant"])
aa_list = list(df["Mutation"])

# amino acids dictionary
aa_dict = {
    "F": ["TTT", "TTC"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "I": ["ATT", "ATC", "ATA"],
    "M": ["ATG"],
    "V": ["GTT", "GTC", "GTA", "GTA"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "Y": ["TAT", "TAC"],
    "H": ["CAT", "CAC"],
    "Q": ["CAA", "CAG"],
    "N": ["AAT", "AAC"],
    "K": ["AAA", "AAG"],
    "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"],
    "C": ["TGT", "TGC"],
    "W": ["TGC"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
    "X": [""],
}


##if X look at position of mutation -> work it out manually -> print anything up to 100 bp
def reverse_complement(input_sequence):
    pairs = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
    }
    reverse_complement = ""
    for index in range(len(input_sequence) - 1, -1, -1):
        base = input_sequence[index]
        complement = pairs[base]
        reverse_complement += complement
    return reverse_complement


def complement(input_sequence):
    pairs = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
    }
    complement = ""
    for base in input_sequence:
        complement += pairs[base]
    return complement


codon_change_dict = {
    8: 4,
    7: 3,
    6: 2,
    5: 1,
    4: 0,
    3: -1,
    2: -2,
    1: -3,
    0: -4,
}


def generate_mutant_sequences(mutation_from, mutation_to, ref_sequences):
    mutant_sequences = []
    mutation_dict = {"A": "A", "G": "G", "C": "C", "T": "T"}
    mutation_dict.update({mutation_from: mutation_to})
    for mutation_position in range(3):
        new_ref_seq = []
        for ref_seq in ref_sequences:
            rs_ = list(ref_seq)
            rs_[mutation_position] = mutation_dict[rs_[mutation_position]]
            new_ref_seq.append("".join(rs_))
        mutant_sequences.append(new_ref_seq)
        # mutant_sequences.append(["".join([mutation_dict[rs] for rs in ref_seq]) for ref_seq in ref_sequences])
        # for ref_seq in ref_sequences:
        #     for rs in ref_seq:
        #         mutation_dict[rs]
        # Generate mutant sequences
        # mutant_sequences.append(
        #     [
        #         seq[:mutation_position] + mutation_to + seq[mutation_position + 1 :]
        #         if seq[mutation_position] == mutation_from
        #         else seq
        #         for seq in ref_sequences
        #     ]
        # )
    mutant_sequences_dict = {}
    print(mutant_sequences)
    print(ref_sequences)
    for i in [0, 1, 2]:
        mutant_sequences_dict[i] = dict(zip(mutant_sequences[i], ref_sequences))
    return mutant_sequences_dict


def match_mutant_sequences(mutant_sequences, target_sequences):
    # Find matching mutant sequences in the target sequences
    matching_indices = [
        index for index, seq in enumerate(target_sequences) if seq in mutant_sequences
    ]
    # Extract matching mutant sequences
    matching_mutant_sequences = [mutant_sequences[index] for index in matching_indices]
    return matching_mutant_sequences


# iterate through list to get sequences of reference and varient
new_seq_final = {}
# for i in range(len(mut_list)):
for i in [3]:
    var = var_list[i]
    ref = ref_list[i]
    ref_to_var_dict = {ref: var}
    if (len(var) == 1) and (len(ref) == 1):
        aa_old, aa_new = aa_list[i].split("/")[0], aa_list[i].split("/")[1]
        codon_old, codon_new = aa_dict[aa_old], aa_dict[aa_new]
        codon_old_rc = [reverse_complement(cod) for cod in codon_old]
        # codon_old_rc = [complement(cod) for cod in codon_old]
        codon_new_rc = [reverse_complement(cod) for cod in codon_new]
        # codon_new_rc = [complement(cod) for cod in codon_new]
        # mutate each codon and see what happens
        codon_old_mutated = generate_mutant_sequences(ref, var, codon_old)
        codon_old_rc_mutated = generate_mutant_sequences(ref, var, codon_old_rc)
        codon_match, codon_rc_match = {}, {}
        for c in codon_old_mutated:
            codon_match[c] = match_mutant_sequences(
                list(codon_old_mutated[c].keys()), codon_new
            )
        for c in codon_old_rc_mutated:
            codon_rc_match[c] = match_mutant_sequences(
                list(codon_old_rc_mutated[c].keys()), codon_new_rc
            )
        mut_region = str(
            reference.get_seq(
                str(chrom_list[0][i]), int(mut_list[i] - 5), int(mut_list[i]) + 5
            )
        )
        mut_region_rc = reverse_complement(mut_region)
        # seq_rc = complement(seq)
        for c in codon_match:
            if len(codon_match[c]) > 0:
                for cod in codon_match[c]:
                    original_codon = codon_old_mutated[c][cod]
                    pos = mut_region.find(original_codon)
                    print(original_codon, pos)
        for c in codon_rc_match:
            if len(codon_rc_match[c]) > 0:
                for cod in codon_rc_match[c]:
                    original_codon = codon_old_rc_mutated[c][cod]
                    pos = mut_region_rc.find(original_codon)
                    print(original_codon, pos)
#         ref_seq_list, mut_seq_list, codon_list = [], [], []
#         ref_seq_rc_list, mut_seq_rc_list, codon_rc_list = [], [], []
#         new_seq = []
#         for x in codon_old:
#             if mut_region.find(str(x)) > 0:
#                 n = mut_region.find(str(x))
#                 nstop = n + 3
#                 for y in codon_new:
#                     x_new = "".join(
#                         [
#                             ref_to_var_dict[xx] if xx in ref_to_var_dict else xx
#                             for xx in x
#                         ]
#                     )
#                     new_mut_region = seq[45:55]
#                     if x_new == y:
#                         new_seq.append(
#                             new_mut_region[:n] + y + new_mut_region[nstop:]
#                         )
#                 if len(new_seq) != 0:
#                     for s in new_seq:
#                         new_seq_long = str(seq)
#                         new_seq_long = new_seq_long[:45] + s + new_seq_long[55:]
#                         new_seq_final.update(
#                             {
#                                 i: new_seq_long[
#                                     4
#                                     + codon_change_dict[n] : 97
#                                     + codon_change_dict[n]
#                                 ]
#                             }
#                         )
#                 else:
#                     for x in codon_old_rc:
#                         if mut_region_rc.find(str(x)) > 0:
#                             n = mut_region_rc.find(str(x))
#                             nstop = n + 3
#                             for y in codon_new_rc:
#                                 x_new = "".join(
#                                     [
#                                         ref_to_var_dict[xx]
#                                         if xx in ref_to_var_dict
#                                         else xx
#                                         for xx in x
#                                     ]
#                                 )
#                                 new_mut_region = seq_rc[45:55]
#                                 if x_new == y:
#                                     new_seq.append(
#                                         new_mut_region[:n]
#                                         + y
#                                         + new_mut_region[nstop:]
#                                     )
#                             for s in new_seq:
#                                 new_seq_long = str(seq_rc)
#                                 new_seq_long = (
#                                     new_seq_long[:45] + s + new_seq_long[55:]
#                                 )
#                                 new_seq_final.update(
#                                     {
#                                         i: new_seq_long[
#                                             4
#                                             + codon_change_dict[n] : 97
#                                             + codon_change_dict[n]
#                                         ]
#                                     }
#                                 )
#         else:
#             new_seq_final.update({i: ""})
#     else:
#         new_seq_final.update({i: ""})
# df2 = pd.DataFrame.from_dict(new_seq_final, orient="index")
# df["new_seq"] = pd.Series(new_seq_final)
# df2.to_csv(
#     "/scratch/project/tcr_neoantigen/results/blood_results_171123/neoantigens/bloodXTumour/Class_I/test.tsv",
#     sep="\t",
#     # index=False,
# )
# if len(l) < 6
# mut =  mut[:45] + var_list[i] + mut[40:]
# seq_list.append(str(seq))
# mut_seq_list.append(mut)
# aa_change = aa_list[i]
# aa = aa_change.partition('/')[2]
# codon = aa_dict[aa]
# if len(codon) > 1:
#     for x in codon:
#         if mut.find(x,37,42) > 0:
#             codon_list.append(x)
# else:
#     codon_list.append(codon)
# print(seq, mut, seq_list, mut_seq_list, aa_change, aa, codon, codon_list)
# break

##Append dataframe of reference and varient sequence
# df['Reference Sequence'] = seq_list
# df['Varient Sequence'] = mut_seq_list
