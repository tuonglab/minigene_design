# From neoepitopes files (class ii) eg. /scratch/project/tcr_neoantigen/results/results_healthyXtumour/neoantigens/sample1/Class_I/sample1_MHC_Class_I_all_epitopes_ccf.tsv
# https://stackoverflow.com/questions/64243197/how-to-query-a-region-in-a-fasta-file-using-biopython
# Get 45bp upstream and downstream of the seq (from neoepitopes) (unique sequence filter)

import pandas as pd
from pyfaidx import Fasta

# import reference and data
ref = Fasta(
    filename="/scratch/project/tcr_neoantigen/resources/references/hg38/gdc/GRCh38.d1.vd1/fasta/GRCh38.d1.vd1.fa"
)
df = pd.read_csv(
    "/scratch/project/tcr_neoantigen/results/blood_results_171123/neoantigens/bloodXTumour/Class_I/bloodXTumour_MHC_Class_I_all_epitopes_ccf_ref_match.tsv",
    sep="\t",
)

# create lists of mutation location and varients
mut_list = list(df["Stop"])
Chrom_list = (list(df["Chromosome"]),)
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


##iterate through list to get sequences of reference and varient
final_list = []
for i in range(len(mut_list)):
    var = var_list[i]
    start_shift = {0: 45, 1: 46, 2: 47}
    end_shift = {0: 47, 1: 46, 2: 45}
    end_shift2 = {0: 47, 1: 48, 2: 49}
    ref_seq_final, mut_seq_final = [], []
    for pos_shift in [0, 1, 2]:
        seq = ref.get_seq(
            str(Chrom_list[0][i]),
            int(mut_list[i] - start_shift[pos_shift]),
            int(mut_list[i]) + end_shift[pos_shift],
        )  # 45 bp upstream and downstream of codon
        ref_seq_final.append(str(seq))
        mut_seq = (
            str(seq[: start_shift[pos_shift]])
            + var
            + str(seq[end_shift2[pos_shift] - 1 :])
        )
        mut_seq_final.append(mut_seq)
    final_list.append(ref_seq_final + mut_seq_final)

dfnew = pd.DataFrame(final_list)
dfnew.columns = [
    "ref_seq_frame_1",
    "ref_seq_frame_2",
    "ref_seq_frame_3",
    "mut_seq_frame_1",
    "mut_seq_frame_2",
    "mut_seq_frame_3",
]

out_df = pd.concat([df, dfnew], axis=1)
out_df.to_csv(
    "/scratch/project/tcr_neoantigen/results/blood_results_171123/neoantigens/bloodXTumour/Class_I/bloodXTumour_MHC_Class_I_all_epitopes_ccf_ref_match_with_93nt.tsv",
    sep="\t",
    index=False,
)
# codon_change_dict = {8:4,
# 7:3,
# 6:2,
# 5:1,
# 4:0,
# 3:-1,
# 2:-2,
# 1:-3,
# 0:-4,}
# seq_list = []
# mut_seq_list = []
# codon_list = []

# refs = ref_list[i]
# ref_to_var_dict = {refs:var}
# aa_change = aa_list[i]
# aa = aa_change.partition('/')[2]
# aa_old = aa_change.split('/')[0]
# aa_new = aa_change.split('/')[1]
# codon_old = aa_dict[aa_old]
# codon_new = aa_dict[aa_new]

# new_seq, new_seq_final = [], []
# for x in codon_old:
#     if mut_region.find(str(x)) > 0:
#         n = mut_region.find(str(x))
#         nstop = n + 3
#         for y in codon_new:
#             x_new = ''.join([ref_to_var_dict[xx] if xx in ref_to_var_dict else xx for xx in x])
#             new_mut_region = str(seq[45:55])
#             if x_new == y:
#                 new_seq.append(new_mut_region[:n] + y + new_mut_region[nstop:])
#         for s in new_seq:
#             new_seq_long = str(seq)
#             new_seq_long = new_seq_long[:45] + s + new_seq_long[55:]
#             new_seq_final.append(new_seq_long[4+codon_change_dict[n]:97+codon_change_dict[n]])
# new_seq_final_final.append(new_seq_final)

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
