file_input_path = "/scratch/project/tcr_neoantigen/results/cSCC_BC_seq_data_10_patients/nextNEOpi/2020239_WO1/analyses/2020239_WO1/05_vep/tables/high_confidence/2020239_WO1_hc_vep.txt"

with open(file_input_path, "r") as file:
    lines = file.readlines()
new_lines = lines[103:]
with open(file_input_path ,"w") as file:
    file.writelines(new_lines)