from pathlib import Path
from pyfaidx import Fasta

from _utils import (
    read_and_filter,
    create_result_list,
    extract_result,
    extract_exon_info,
    filter_exon_pos,
)

# define head paths
HG38FOLDER = Path("/scratch/project/tcr_neoantigen/resources/references/hg38")
INPUTFOLDER = Path(
    "/scratch/project/tcr_neoantigen/results/cSCC_BC_seq_data_10_patients/nextNEOpi"
)

# import reference and data
fasta_file = HG38FOLDER / "gdc" / "GRCh38.d1.vd1" / "fasta" / "GRCh38.d1.vd1.fa"
gtf_file = HG38FOLDER / "annotation" / "gencode.v33.primary_assembly.annotation.gtf"
refgen = Fasta(filename=fasta_file)

samples = [
    # "2020135",
    # "2020239_WO1",
    # "2020246_NO1",
    # "2020260_WO1",
    # "2020281_WO1",
    # "2021111_MO1",
    "DES001",
    # "DES002",
    # "DES002_001",
    # "DES002_002",
    # "DES010",
]

results = {}
for sample in samples:
    file_input_path = (
        INPUTFOLDER
        / sample
        / "analyses"
        / sample
        / "05_vep"
        / "tables"
        / "high_confidence"
        / f"{sample}_hc_vep.txt"
    )
    results[sample] = read_and_filter(file_input_path)

for sample in samples:
    print(sample, results[sample].shape)

for sample in samples:
    print(sample, results[sample].VARIANT_CLASS.unique())

# let's setup with DES001 as the case as it contains all classes we want
df = results["DES001"].copy()

test_dict = create_result_list(df)
# from _utils import find
# >>> find(test_dict["variant_class"], "insertion")
# [2691, 2719, 3920]
# >>> find(test_dict["variant_class"], "deletion")
# [463, 661, 926, 1060, 1242, 1958, 2321, 2330, 2331, 2332, 2566, 2684, 3107, 3194, 3259]
# >>> find(test_dict["variant_class"], "substitution")
# [6, 9, 12, 14, 28, 31, 34, 37, 39, 48, 51, 63, 68, 80, 87, 88, 91, 94, 95, 108, 113, 136, 171, 183, 191, 201, 204, 211, 220, 222, 230, 231, 233, 255, 258, 263, 270, 280, 281, 288, 293, 306, 316, 328, 329, 352, 354, 376, 387, 388, 391, 398, 404, 407, 435, 446, 470, 476, 485, 486, 488, 496, 501, 505, 507, 518, 520, 535, 541, 543, 545, 569, 572, 579, 585, 604, 626, 646, 656, 671, 705, 718, 738, 739, 745, 774, 778, 783, 791, 798, 806, 807, 808, 816, 824, 846, 850, 867, 872, 875, 909, 915, 925, 933, 935, 946, 960, 962, 972, 977, 1011, 1055, 1061, 1066, 1069, 1091, 1101, 1103, 1110, 1127, 1131, 1194, 1232, 1244, 1251, 1256, 1262, 1271, 1294, 1295, 1305, 1315, 1316, 1321, 1323, 1326, 1332, 1347, 1356, 1377, 1401, 1405, 1406, 1409, 1428, 1449, 1470, 1473, 1481, 1501, 1517, 1525, 1526, 1527, 1534, 1545, 1609, 1622, 1623, 1647, 1655, 1658, 1687, 1693, 1701, 1706, 1728, 1782, 1790, 1813, 1834, 1837, 1841, 1843, 1847, 1856, 1902, 1918, 1924, 1928, 1929, 1935, 1939, 1940, 1943, 1977, 1989, 1992, 1994, 1998, 2007, 2009, 2011, 2022, 2027, 2053, 2057, 2060, 2079, 2102, 2106, 2119, 2121, 2157, 2159, 2165, 2175, 2186, 2187, 2191, 2219, 2224, 2229, 2231, 2257, 2262, 2267, 2278, 2279, 2300, 2302, 2315, 2316, 2318, 2328, 2339, 2373, 2376, 2385, 2408, 2419, 2426, 2428, 2476, 2486, 2510, 2511, 2516, 2526, 2535, 2537, 2557, 2559, 2583, 2589, 2634, 2636, 2643, 2647, 2649, 2663, 2664, 2674, 2726, 2779, 2805, 2807, 2827, 2856, 2862, 2868, 2880, 2882, 2883, 2884, 2886, 2887, 2895, 2896, 2898, 2901, 2903, 2905, 2906, 2934, 2949, 2978, 2990, 2995, 3000, 3020, 3032, 3039, 3040, 3058, 3066, 3068, 3069, 3077, 3080, 3094, 3096, 3108, 3114, 3120, 3158, 3164, 3171, 3181, 3218, 3233, 3238, 3256, 3270, 3287, 3296, 3331, 3334, 3342, 3343, 3349, 3352, 3356, 3360, 3369, 3371, 3377, 3384, 3387, 3400, 3401, 3403, 3407, 3410, 3417, 3419, 3423, 3426, 3427, 3430, 3432, 3436, 3449, 3459, 3465, 3473, 3474, 3487, 3512, 3522, 3524, 3528, 3529, 3532, 3533, 3535, 3538, 3547, 3549, 3551, 3556, 3568, 3577, 3584, 3586, 3591, 3596, 3597, 3611, 3617, 3620, 3623, 3628, 3630, 3632, 3634, 3635, 3650, 3651, 3654, 3658, 3661, 3664, 3671, 3680, 3682, 3683, 3691, 3695, 3710, 3734, 3736, 3740, 3745, 3746, 3775, 3790, 3800, 3829, 3833, 3836, 3862, 3866, 3883, 3885, 3887, 3902, 3903, 3914, 3927, 3937, 3941, 3944, 3948, 3967, 3974, 3975, 3985, 3988, 3994, 3995, 4009, 4010, 4011, 4017, 4019, 4020, 4023, 4032, 4057, 4092, 4108, 4112, 4116, 4129, 4130, 4131, 4134, 4139, 4152, 4153, 4156]

test = extract_result(test_dict, 3920)

ref, var, aas, chrom, mut, var_class, gene, gene_name, strand = (
    test["reference"],
    test["variant"],
    test["amino_acid"],
    test["chromosome"],
    test["mutation_location"],
    test["variant_class"],
    test["gene_id"],
    test["gene_symbol"],
    test["strand"],
)

exon_infos = extract_exon_info(gtf_file)
if var_class in ["deletion", "insertion", "substitution"]:
    # we need to create 45 bp overlapping sequences across all exon sequences
    # for example, if we have a gene like so with a framsehifting indel (or stop codon) somewhere in the exon (marked by --); big N are exon regions and small n are inton regions:
    # NNNNNNNNNNNNNNNNNNNNNnnnnnnnNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNnnnnnnnnnnnnnnnNNNNNNNNNNNN--NNNNNNNNNNNNNNnnnnnnnnnnnnnnnnnnnnnnNNNNNNNNNNNNNNNNNNNNNNNnnnnnnnnnnnnnnnNNNNNNNNNNNNNNNNNNNNNNNNNNNNnnnnnnnnnnnnnnnNNNNNNNNNNNNNNNNNNNNNNNNNNnnnnnnnnnnnnnnnnnnnnnnNNNNNNNNNNNNNNNNNNNNNNNnnnnnNNNNNNNNNNNNNNNNNNNNnnnnnNNNNNNNNNNNNNNNNNNNNnnnnnnNNNNNNNNNNNNNNNNNNNNnnnnn
    # so if we get rid of the introns, the sequence should look like this (i'm just using | to show where the boundaries are):
    # NNNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNN--NNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNN
    # Next, we need to consider the following:
    # 1) find the sequences flanking the indel/stop so that the indel/stop is in the middle.
    # 2) We should only use exonic sequences, so we need to obtain the exon positions. the earlier filtering would have trim the variants to only contain such mutations.
    # 3) we would need to create the sequences, by overlaying the last 45bp:
    # transcript:   NNNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNN--NNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNN
    # 1st seq:                        NNN|NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNN--NNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNNNNN|NNNNNNNNN
    # 2nd seq:                                                                          NNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNNNNNNNN|NNN
    # 3rd seq:                                                                                                                            NNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNN|NNNNNNNN
    # 4th seq:                                                                                                                                                             NNNNNNNNNN|NNNNNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNN|NNNNNNNNNNNNNNNNNNNN
    # for the last sequence, we extend leftwards instead of right
    # If possible, we should try to ensure that open reading frame before extending the sequences.
    # We can do this by checking the codon. If not possible, we will have to create 3 separate frames.
    # if it's a deletion, then we don't need to care too much and just extend until the sequence is 93 nt
    # if it's insertion, we extend left and right slowly until the length is 93nt..?
    # if it's stop, the sequence must be exactly mosition 45 before, 46, 47, 48 and 45 after.
    # if the mutation is really on the left hand (5' side), then the mutation don't need to be placed in the middle and we just extend till it's 93nt
    # if there is a stop, Jaz also wants to find the next AUG sequence and see if a new protein is made. We will worry about it when we reach that stage

    # get only relevant exons
    filtered_exon_pos = filter_exon_pos(exon_infos[gene], mut)
    # if we examined the first window of the filtered exon positions, we check if the start position for the mutation and the exon can be extended by 45bp.
    # we also check if the end position is wihin the exon. hopefully it always is.
    # as there may be multiple protein isoforms, we have to iterate for each
    for prot_id in filtered_exon_pos:
        exon_start_diff = mut - filtered_exon_pos[prot_id][0][0]
        if exon_start_diff >= 45:
            if var_class == "insertion":
                # left side
                start_seq_loc = mut - 45
                left_flank = refgen.get_seq(
                    chrom,
                    start_seq_loc,
                    mut - 1,
                )
                # middle
                middle_ref = "" if var_class == "insertion" else ref
                middle_var = "" if var_class == "deletion" else var
                # right side
                end_seq_loc_ref = (
                    mut + 47
                )  # because the location of the mutation counts as 1
                right_flank_ref = refgen.get_seq(chrom, mut, end_seq_loc_ref)
                end_seq_loc_var = mut + (47 - len(var))  # this is the pad
                exon_end_diff = filtered_exon_pos[prot_id][0][1] - end_seq_loc_var
                if exon_end_diff > 0:
                    right_flank_var = refgen.get_seq(chrom, mut, end_seq_loc_var)
                else:
                    right_flank_var = refgen.get_seq(
                        chrom, mut, end_seq_loc + exon_end_diff
                    )
                # get right hang
                right_hang_ref = str(
                    refgen.get_seq(
                        chrom, end_seq_loc_ref + 1, filtered_exon_pos[prot_id][0][1]
                    )
                )
                right_hang_var = str(
                    refgen.get_seq(
                        chrom, end_seq_loc_var + 1, filtered_exon_pos[prot_id][0][1]
                    )
                )
                for n in range(1, len(filtered_exon_pos[prot_id])):
                    right_hang_ref += str(
                        refgen.get_seq(
                            chrom,
                            filtered_exon_pos[prot_id][n][0],
                            filtered_exon_pos[prot_id][n][1],
                        )
                    )
                    right_hang_var += str(
                        refgen.get_seq(
                            chrom,
                            filtered_exon_pos[prot_id][n][0],
                            filtered_exon_pos[prot_id][n][1],
                        )
                    )
                # split the right hangs into equal 93 nt, with the last one just extending from the left

# else:
#     # current workflow
#     pass


# aa_old, aa_new = aas.split("/")[0], aas.split("/")[1]
# start_shift = {0: 45, 1: 46, 2: 47}
# end_shift = {0: 47, 1: 46, 2: 45}
# end_shift2 = {0: 47, 1: 48, 2: 49}
# ref_seq_final, mut_seq_final = [], []
# ref_seq_rc_final, mut_seq_rc_final = [], []
# ref_seq_final_aa, mut_seq_final_aa = [], []
# ref_seq_rc_final_aa, mut_seq_rc_final_aa = [], []
# match_final = []
# for pos_shift in [0, 1, 2]:
#     if var == "-":
#         # its a deletion
#         ref_seq = str(
#             refgen.get_seq(
#                 chrom,
#                 locs[0] - 45 + pos_shift + 1,
#                 locs[1] + 45 + pos_shift,  # if deletion
#             )
#         )
#         print(ref_seq)  # 45 bp upstream and downstream of codon
#         # so the windows are
#         # (48) 45 + 3 + 45 + 0 (45)
#         # (47) 44 + 3 + 45 + 1 (46)
#         # (46) 43 + 3 + 45 + 2 (47)
#         mut_seq = (
#             ref_seq[: start_shift[pos_shift]] + ref_seq[end_shift2[pos_shift] - 1 :]
#         )
#         print(mut_seq)


# ref_seq = str(
#     refgen.get_seq(
#         chrom,
#         exon_pos[0],
#         locs[1] + 45,  # if deletion
#     )
# )
