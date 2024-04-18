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
# drop to unique Stop, Reference, Variant rows
df = df.drop_duplicates(["Stop", "Reference", "Variant"])
df.to_csv(
    "/scratch/project/tcr_neoantigen/results/blood_results_171123/neoantigens/bloodXTumour/Class_I/bloodXTumour_MHC_Class_I_all_epitopes_ccf_ref_match_small.tsv",
    sep="\t",
    index=False,
)
codon_dict = {
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "C": ["TGT", "TGC"],
    "W": ["TGG"],
    "E": ["GAA", "GAG"],
    "D": ["GAT", "GAC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "N": ["AAT", "AAC"],
    "M": ["ATG"],
    "K": ["AAA", "AAG"],
    "Y": ["TAT", "TAC"],
    "I": ["ATT", "ATC", "ATA"],
    "Q": ["CAA", "CAG"],
    "F": ["TTT", "TTC"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "*": ["TAA", "TAG", "TGA"],
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
    "H": ["CAT", "CAC"],
}

mut_list = list(df["Stop"])
chrom_list = (list(df["Chromosome"]),)
ref_list = list(df["Reference"])
var_list = list(df["Variant"])
aa_change_list = list(df["Mutation"])


def reverse_complement(sequence):
    complement_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complement_dict[base] for base in reversed(sequence))


seq_possible_result, mut_possible_result = [], []
for i in range(len(mut_list)):
    var = var_list[i]
    aa_change = aa_change_list[i]
    aa_ref = aa_change.split("/")[0]
    aa_mut = aa_change.split("/")[1]
    seq_possible, mut_possible = [], []
    if (len(aa_ref) == 1) & (len(aa_mut) == 1):  # focus on SNV for now
        if aa_mut != "X":  # unknown change?
            codon_ref = codon_dict[aa_ref]
            codon_mut = codon_dict[aa_mut]
            # get the 3 possible codons:
            # x - -
            # - X -
            # - - X
            start_pos_shift_dict = {0: 0, 1: -1, 2: -2}
            start_pos_shift_dict_for_mut = {0: 0, 1: 1, 2: 2}
            end_pos_shift_dict = {0: 2, 1: 1, 2: 0}
            end_pos_shift_dict_for_mut = {0: 1, 1: 2, 2: 3}
            for pos_shift in [0, 1, 2]:
                seq = str(
                    ref.get_seq(
                        str(chrom_list[0][i]),
                        mut_list[i] + start_pos_shift_dict[pos_shift],
                        mut_list[i] + end_pos_shift_dict[pos_shift],
                    )
                )
                mut = (
                    seq[: start_pos_shift_dict_for_mut[pos_shift]]
                    + var
                    + seq[end_pos_shift_dict_for_mut[pos_shift] :]
                )
                seq_possible.append(seq in codon_ref)
                mut_possible.append(mut in codon_mut)
                rev_seq = reverse_complement(seq)
                rev_mut = reverse_complement(mut)
                seq_possible.append(rev_seq in codon_ref)
                mut_possible.append(rev_mut in codon_mut)
    seq_possible_result.append(1 if any(seq_possible) else 0)
    mut_possible_result.append(1 if any(mut_possible) else 0)

possible = pd.DataFrame([seq_possible_result, mut_possible_result]).T
possible.columns = ["ref_possible", "mut_possible"]
possible.to_csv(
    "/scratch/project/tcr_neoantigen/results/blood_results_171123/neoantigens/bloodXTumour/Class_I/bloodXTumour_MHC_Class_I_check_possible.tsv",
    sep="\t",
    index=False,
)
