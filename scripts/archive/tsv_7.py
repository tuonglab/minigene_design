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
df = df.drop_duplicates(["Stop", "Reference", "Variant"])
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


def translate(seq):
    nnn_table = {
        "TTT": "F",
        "TCT": "S",
        "TAT": "Y",
        "TGT": "C",
        "TTC": "F",
        "TCC": "S",
        "TAC": "Y",
        "TGC": "C",
        "TTA": "L",
        "TCA": "S",
        "TAA": "*",
        "TGA": "*",
        "TTG": "L",
        "TCG": "S",
        "TAG": "*",
        "TGG": "W",
        "CTT": "L",
        "CCT": "P",
        "CAT": "H",
        "CGT": "R",
        "CTC": "L",
        "CCC": "P",
        "CAC": "H",
        "CGC": "R",
        "CTA": "L",
        "CCA": "P",
        "CAA": "Q",
        "CGA": "R",
        "CTG": "L",
        "CCG": "P",
        "CAG": "Q",
        "CGG": "R",
        "ATT": "I",
        "ACT": "T",
        "AAT": "N",
        "AGT": "S",
        "ATC": "I",
        "ACC": "T",
        "AAC": "N",
        "AGC": "S",
        "ATA": "I",
        "ACA": "T",
        "AAA": "K",
        "AGA": "R",
        "ATG": "M",
        "ACG": "T",
        "AAG": "K",
        "AGG": "R",
        "GTT": "V",
        "GCT": "A",
        "GAT": "D",
        "GGT": "G",
        "GTC": "V",
        "GCC": "A",
        "GAC": "D",
        "GGC": "G",
        "GTA": "V",
        "GCA": "A",
        "GAA": "E",
        "GGA": "G",
        "GTG": "V",
        "GCG": "A",
        "GAG": "E",
        "GGG": "G",
    }
    # two loops, outer one to loop over the list of string sequences
    # inner one loops over each sequence
    nnn_aa_seq = []
    codons = [seq[i : i + 3] for i in range(0, len(seq), 3)]
    for c in codons:
        nnn_aa_seq.append(nnn_table[c])
    return "".join(nnn_aa_seq)


def match15(aa: str, match: str):
    """This amino acid position is the position where we are expecting the mutation!"""
    return aa[15] == match


##iterate through list to get sequences of reference and varient
final_list = []
for i in range(len(mut_list)):
    var = var_list[i]
    aa_old, aa_new = aa_list[i].split("/")[0], aa_list[i].split("/")[1]
    start_shift = {0: 45, 1: 46, 2: 47}
    end_shift = {0: 47, 1: 46, 2: 45}
    end_shift2 = {0: 47, 1: 48, 2: 49}
    ref_seq_final, mut_seq_final = [], []
    ref_seq_rc_final, mut_seq_rc_final = [], []
    ref_seq_final_aa, mut_seq_final_aa = [], []
    ref_seq_rc_final_aa, mut_seq_rc_final_aa = [], []
    match_final = []
    for pos_shift in [0, 1, 2]:
        ref_seq = str(
            ref.get_seq(
                str(Chrom_list[0][i]),
                int(mut_list[i] - start_shift[pos_shift]),
                int(mut_list[i]) + end_shift[pos_shift],
            )
        )  # 45 bp upstream and downstream of codon
        mut_seq = (
            ref_seq[: start_shift[pos_shift]]
            + var
            + ref_seq[end_shift2[pos_shift] - 1 :]
        )
        # reverse complement the whole sequence
        ref_seq_rc = reverse_complement(ref_seq)
        mut_seq_rc = reverse_complement(mut_seq)
        # translate to aa
        ref_seq_aa = translate(ref_seq)
        mut_seq_aa = translate(mut_seq)
        ref_seq_rc_aa = translate(ref_seq_rc)
        mut_seq_rc_aa = translate(mut_seq_rc)
        # query the original and new aa to see if they match
        ref_seq_aa_match = match15(ref_seq_aa, aa_old)
        mut_seq_aa_match = match15(mut_seq_aa, aa_new)
        ref_seq_rc_aa_match = match15(ref_seq_rc_aa, aa_old)
        mut_seq_rc_aa_match = match15(mut_seq_rc_aa, aa_new)
        # append everything.
        if not ref_seq_aa_match:
            ref_seq, ref_seq_aa, mut_seq, mut_seq_aa = "", "", "", ""
        else:
            if not mut_seq_aa_match:
                ref_seq, ref_seq_aa, mut_seq, mut_seq_aa = "", "", "", ""
        if not ref_seq_rc_aa_match:
            ref_seq_rc, ref_seq_rc_aa, mut_seq_rc, mut_seq_rc_aa = "", "", "", ""
        else:
            if not mut_seq_rc_aa_match:
                ref_seq_rc, ref_seq_rc_aa, mut_seq_rc, mut_seq_rc_aa = "", "", "", ""
        ref_seq_final.append(ref_seq) if ref_seq != "" else None
        mut_seq_final.append(mut_seq) if mut_seq != "" else None
        ref_seq_rc_final.append(ref_seq_rc) if ref_seq_rc != "" else None
        mut_seq_rc_final.append(mut_seq_rc) if mut_seq_rc != "" else None
        ref_seq_final_aa.append(ref_seq_aa) if ref_seq_aa != "" else None
        mut_seq_final_aa.append(mut_seq_aa) if mut_seq_aa != "" else None
        ref_seq_rc_final_aa.append(ref_seq_rc_aa) if ref_seq_rc_aa != "" else None
        mut_seq_rc_final_aa.append(mut_seq_rc_aa) if mut_seq_rc_aa != "" else None
        if ref_seq_aa_match and mut_seq_aa_match:
            match_final.append("+")
        if ref_seq_rc_aa_match and mut_seq_rc_aa_match:
            match_final.append("-")
    final_list.append(
        [aa_old]
        + [aa_new]
        + ref_seq_final
        + mut_seq_final
        + ref_seq_rc_final
        + mut_seq_rc_final
        + ref_seq_final_aa
        + mut_seq_final_aa
        + ref_seq_rc_final_aa
        + mut_seq_rc_final_aa
        + match_final
    )

dfnew = pd.DataFrame(final_list)
dfnew.columns = [
    "ref",
    "alt",
    "ref_seq_frame_1",
    "ref_seq_frame_2",
    "ref_seq_frame_3",
    "mut_seq_frame_1",
    "mut_seq_frame_2",
    "mut_seq_frame_3",
    "ref_seq_reverse_complement_frame_1",
    "ref_seq_reverse_complement_frame_2",
    "ref_seq_reverse_complement_frame_3",
    "mut_seq_reverse_complement_frame_1",
    "mut_seq_reverse_complement_frame_2",
    "mut_seq_reverse_complement_frame_3",
    "ref_seq_aa_frame_1",
    "ref_seq_aa_frame_2",
    "ref_seq_aa_frame_3",
    "mut_seq_aa_frame_1",
    "mut_seq_aa_frame_2",
    "mut_seq_aa_frame_3",
    "ref_seq_reverse_complement_aa_frame_1",
    "ref_seq_reverse_complement_aa_frame_2",
    "ref_seq_reverse_complement_aa_frame_3",
    "mut_seq_reverse_complement_aa_frame_1",
    "mut_seq_reverse_complement_aa_frame_2",
    "mut_seq_reverse_complement_aa_frame_3",
    "ref_seq_aa_frame_1_match",
    "ref_seq_aa_frame_2_match",
    "ref_seq_aa_frame_3_match",
    "mut_seq_aa_frame_1_match",
    "mut_seq_aa_frame_2_match",
    "mut_seq_aa_frame_3_match",
    "ref_seq_reverse_complement_aa_frame_1_match",
    "ref_seq_reverse_complement_aa_frame_2_match",
    "ref_seq_reverse_complement_aa_frame_3_match",
    "mut_seq_reverse_complement_aa_frame_1_match",
    "mut_seq_reverse_complement_aa_frame_2_match",
    "mut_seq_reverse_complement_aa_frame_3_match",
]

dfnew = dfnew[
    [
        "ref",
        "alt",
        "ref_seq_frame_1",
        "mut_seq_frame_1",
        "ref_seq_frame_2",
        "mut_seq_frame_2",
        "ref_seq_frame_3",
        "mut_seq_frame_3",
        "ref_seq_reverse_complement_frame_1",
        "mut_seq_reverse_complement_frame_1",
        "ref_seq_reverse_complement_frame_2",
        "mut_seq_reverse_complement_frame_2",
        "ref_seq_reverse_complement_frame_3",
        "mut_seq_reverse_complement_frame_3",
        "ref_seq_aa_frame_1",
        "mut_seq_aa_frame_1",
        "ref_seq_aa_frame_2",
        "mut_seq_aa_frame_2",
        "ref_seq_aa_frame_3",
        "mut_seq_aa_frame_3",
        "ref_seq_reverse_complement_aa_frame_1",
        "mut_seq_reverse_complement_aa_frame_1",
        "ref_seq_reverse_complement_aa_frame_2",
        "mut_seq_reverse_complement_aa_frame_2",
        "ref_seq_reverse_complement_aa_frame_3",
        "mut_seq_reverse_complement_aa_frame_3",
        "ref_seq_aa_frame_1_match",
        "mut_seq_aa_frame_1_match",
        "ref_seq_aa_frame_2_match",
        "mut_seq_aa_frame_2_match",
        "ref_seq_aa_frame_3_match",
        "mut_seq_aa_frame_3_match",
        "ref_seq_reverse_complement_aa_frame_1_match",
        "mut_seq_reverse_complement_aa_frame_1_match",
        "ref_seq_reverse_complement_aa_frame_2_match",
        "mut_seq_reverse_complement_aa_frame_2_match",
        "ref_seq_reverse_complement_aa_frame_3_match",
        "mut_seq_reverse_complement_aa_frame_3_match",
    ]
]

# out_df = pd.concat([df, dfnew], axis=1)
dfnew.to_csv(
    "/scratch/project/tcr_neoantigen/results/blood_results_171123/neoantigens/bloodXTumour/Class_I/bloodXTumour_MHC_Class_I_all_epitopes_ccf_ref_match_with_93nt.tsv",
    sep="\t",
    # index=False,
)
