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


##append to df dont need to drop duplicates
# drop duplicates
df.drop_duplicates(subset=["Chromosome", "Start", "Stop"], keep="first", inplace=True)

# create lists of mutation location and varients
mut_list = list(df["Start"])
Chrom_list = (list(df["Chromosome"]),)
var_list = list(df["Variant"])


##iterate through list to get sequences of reference and varient
location_list = []
seq_list = []
mut_seq_list = []
for i in range(len(mut_list)):
    seq = ref.get_seq(
        str(Chrom_list[0][i]), int(mut_list[i] - 45), int(mut_list[i]) + 45
    )  # 45 bp upstream and downstream
    mut = str(seq)
    mut = mut[:39] + var_list[i] + mut[40:]
    location = [str(Chrom_list[0][i]), int(mut_list[i] - 45), int(mut_list[i]) + 45]
    seq_list.append(str(seq))
    mut_seq_list.append(mut)
    location_list.append(location)


##create dataframe of location, reference and varient sequence
data = {
    "Location": location_list,
    "Reference Sequence": seq_list,
    "Varient Sequence": mut_seq_list,
}
df2 = pd.DataFrame(data, columns=["Location", "Reference Sequence", "Varient Sequence"])
print(df2)
