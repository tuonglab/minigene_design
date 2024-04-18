import re
import shutil

import pandas as pd

from collections import defaultdict
from pathlib import Path

FINAL_SEQ_LEN = 93


def read_and_filter(
    file_input_path: str | Path, keep_from_lines: int = 103
) -> pd.DataFrame:
    """
    _summary_

    Parameters
    ----------
    file_input_path : str | Path

    keep_from_lines : int, optional
        _description_, by default 103

    Returns
    -------
    pd.DataFrame
        _description_
    """
    file_input_path = (
        Path(file_input_path)
        if not isinstance(file_input_path, Path)
        else file_input_path
    )
    try:
        df = pd.read_csv(file_input_path, sep="\t")
    except:  # from Alicia's clean_up_files.py
        shutil.copy(
            file_input_path,
            file_input_path.parent
            / (file_input_path.stem + "_old" + file_input_path.suffix),
        )
        with open(file_input_path, "r") as file:
            lines = file.readlines()
        new_lines = lines[keep_from_lines:]
        with open(file_input_path, "w") as file:
            file.writelines(new_lines)
        df = pd.read_csv(file_input_path, sep="\t")
    df = df.drop_duplicates(["Location"])
    df["Reference"] = [i.split("_")[2].split("/")[0] for i in df["#Uploaded_variation"]]
    df = df[df["IMPACT"].isin(["LOW", "MODERATE", "HIGH"])]
    # minimal filter
    df["Keep"] = False
    for i, r in df.iterrows():
        if re.search("/", r["Amino_acids"]):
            df.at[i, "Keep"] = True
        if len(r["Reference"]) > 1:
            df.at[i, "Keep"] = True
        if r["SYMBOL"] == "-":
            df.at[i, "Keep"] = False
    df = df[df.Keep]
    df["chromosome"] = [r.split("_")[0] for r in df["#Uploaded_variation"]]
    df["reference"] = [r.split("_")[2].split("/")[0] for r in df["#Uploaded_variation"]]
    df["variant"] = [r.split("_")[2].split("/")[1] for r in df["#Uploaded_variation"]]
    df["aa_ref"] = [r.split("/")[0] for r in df["Amino_acids"]]
    df["aa_var"] = [
        r.split("/")[1] if re.search("/", r) else r for r in df["Amino_acids"]
    ]
    df["codon_ref"] = [r.split("/")[0] for r in df["Codons"]]
    df["codon_var"] = [
        r.split("/")[1] if re.search("/", r) else r for r in df["Codons"]
    ]
    df["mutation_location"] = [
        (
            int(r.split(":")[1].split("-")[0])
            if re.search("-", r)
            else int(r.split(":")[1])
        )
        for r in df["Location"]
    ]
    df["mutation_location2"] = [
        (
            int(r.split(":")[1].split("-")[1])
            if re.search("-", r)
            else int(r.split(":")[1])
        )
        for r in df["Location"]
    ]
    for i, r in df.iterrows():
        if r["aa_ref"] == r["aa_var"]:
            df.at[i, "Keep"] = False
    df = df[df.Keep]
    return df


def complementary_sequence(dna_sequence):
    complement_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    complementary_seq = "".join(
        complement_dict.get(base, base) for base in dna_sequence
    )
    return complementary_seq


def create_result_list(df: pd.DataFrame) -> list:
    chrom_list = list(df["chromosome"])
    ref_list = list(df["reference"])
    var_list = list(df["variant"])
    mut_list = list(df["mutation_location"])
    mut_list2 = list(df["mutation_location2"])
    aa_ref_list = list(df["aa_ref"])
    aa_var_list = list(df["aa_var"])
    codon_ref_list = list(df["codon_ref"])
    codon_var_list = list(df["codon_var"])
    var_class_list = list(df["VARIANT_CLASS"])
    gene_list = list(df["Gene"])
    ensp_list = list(df["ENSP"])
    enst_list = list(df["Feature"])
    gene_name_list = list(df["SYMBOL"])
    strand_list = list(df["STRAND"])
    return {
        "reference": ref_list,
        "variant": var_list,
        "amino_acid_ref": aa_ref_list,
        "amino_acid_var": aa_var_list,
        "codon_ref": codon_ref_list,
        "codon_var": codon_var_list,
        "chromosome": chrom_list,
        "mutation_location": mut_list,
        "mutation_location2": mut_list2,
        "variant_class": var_class_list,
        "gene_id": gene_list,
        "protein_id": ensp_list,
        "transcript_id": enst_list,
        "gene_symbol": gene_name_list,
        "strand": strand_list,
    }


def extract_result(df_dict: dict, index: int) -> dict:
    return {
        "reference": df_dict["reference"][index],
        "variant": df_dict["variant"][index],
        "amino_acid_ref": df_dict["amino_acid_ref"][index],
        "amino_acid_var": df_dict["amino_acid_var"][index],
        "codon_ref": df_dict["codon_ref"][index],
        "codon_var": df_dict["codon_var"][index],
        "chromosome": df_dict["chromosome"][index],
        "mutation_location": df_dict["mutation_location"][index],
        "mutation_location2": df_dict["mutation_location2"][index],
        "variant_class": df_dict["variant_class"][index],
        "gene_id": df_dict["gene_id"][index],
        "protein_id": df_dict["protein_id"][index],
        "transcript_id": df_dict["transcript_id"][index],
        "gene_symbol": df_dict["gene_symbol"][index],
        "strand": df_dict["strand"][index],
    }


def extract_exon_info(gtf_file):
    # Parse the GTF file to get exon positions
    gene_info = get_exon_info(gtf_file)
    # Extract exon sequences
    exon_positions = defaultdict(lambda: defaultdict(list))
    for gene_id, (_, chrom, positions) in gene_info.items():
        if chrom.startswith("chr"):
            new_gene_id = gene_id.split(".")[0]
            for protein_id, exon_pos_list in positions.items():
                if protein_id is not None:
                    exon_positions[new_gene_id][protein_id].extend(exon_pos_list)
    return exon_positions


def get_exon_info(gtf_file):
    exon_positions = {}
    gene_info = {}
    # Read the GTF file and extract exon positions
    with open(gtf_file, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            feature_type = fields[2]
            chrom = fields[0]
            if feature_type == "exon":
                gene_id, gene_name, gene_type, protein_id = get_gene_info(fields[8])
                if gene_type == "protein_coding":
                    start, end = map(int, fields[3:5])
                    if gene_id not in exon_positions:
                        exon_positions[gene_id] = defaultdict(list)
                    exon_positions[gene_id][protein_id].append((start, end))
                    gene_info[gene_id] = (
                        gene_name,
                        chrom,
                        exon_positions[gene_id],
                    )
    return gene_info


def get_gene_info(attributes):
    # Parse the attributes field in the GTF file to get the gene_id and gene_name
    attributes = attributes.split(";")
    gene_id, gene_name, gene_type, protein_id = None, None, None, None
    for attribute in attributes:
        if attribute != "":
            key, value = attribute.strip().split(" ")
            if key == "gene_type":
                gene_type = value.strip('"')
            if key == "gene_id":
                gene_id = value.strip('"')
            if key == "gene_name":
                gene_name = value.strip('"')
            if key == "protein_id":
                protein_id = value.strip('"')
    return gene_id, gene_name, gene_type, protein_id


def closest_smaller_number(exon_pos, loc):
    closest_tuple = None
    for exon in exon_pos:
        if loc > exon[0] and loc < exon[1]:
            closest_tuple = exon
            break
    return closest_tuple


def closest_larger_number(exon_pos, loc):
    closest_tuple = None
    for exon in exon_pos:
        if exon[0] < loc and not exon[1] < loc:
            closest_tuple = exon
            break
    return closest_tuple


def filter_exon_pos(exon_pos, loc, mut_info, strand):
    result = {}
    for prot_id in exon_pos:
        if prot_id.split(".")[0] == mut_info["protein_id"]:
            if strand == "1":
                closest = closest_smaller_number(exon_pos[prot_id], loc)
            else:
                closest = closest_larger_number(exon_pos[prot_id], loc)
            if closest is not None:
                result[prot_id] = exon_pos[prot_id][exon_pos[prot_id].index(closest) :]
    return result


def find(lst, search):
    return [i for i, x in enumerate(lst) if x == search]


def generate_windows(sequence, window_size=93, overlap=45):
    windows = []
    start = 0
    end = window_size
    while start < len(sequence):
        window = sequence[start:end]
        if len(window) == window_size:
            windows.append(window)
            start += window_size - overlap
            end = min(start + window_size, len(sequence))
        else:
            break
    last_window_size = min(window_size, len(sequence) - start)
    if last_window_size % 3 != 0 or last_window_size < window_size:
        right_trim = last_window_size % 3  # Calculate how much to trim
        last_window_size -= right_trim  # Trim the right side
        left_extension = window_size - last_window_size
        start -= left_extension
    end = start + window_size  # Adjust end
    last_window = sequence[start:end]
    windows.append(last_window)
    return windows


def print_windows(sequence, windows, window_size=93, overlap=45):
    for i, window in enumerate(windows[:-1]):
        left_space = " " * (i * (window_size - overlap))
        right_space = " " * (window_size - overlap - len(window))
        print(left_space + window + right_space)

    # Print the last window with the appropriate spacing
    last_left_space = " " * (len(sequence) - window_size)
    print(last_left_space + windows[-1])


def get_right_hang(
    first_seq,
    chrom,
    prot_id,
    first_seq_end_loc,
    filtered_exon_pos,
    fasta,
):
    right_hang = first_seq[-45:]
    right_hang_pad = fasta.get_seq(
        chrom, first_seq_end_loc + 1, filtered_exon_pos[prot_id][0][1]
    )
    right_hang += str(right_hang_pad)
    for n in range(1, len(filtered_exon_pos[prot_id])):
        add_hang = fasta.get_seq(
            chrom,
            filtered_exon_pos[prot_id][n][0],
            filtered_exon_pos[prot_id][n][1],
        )
        right_hang += str(add_hang)
    # generate the overlapping windows
    overlap = generate_windows(right_hang)
    return overlap


def get_left_hang(
    first_seq,
    chrom,
    prot_id,
    first_seq_start_loc,
    filtered_exon_pos,
    fasta,
):
    left_hang = first_seq[:45]  # Get the first 45 characters from the start
    left_hang_pad = fasta.get_seq(
        chrom, filtered_exon_pos[prot_id][0][0], first_seq_start_loc - 1
    )
    left_hang = str(left_hang_pad) + left_hang
    for n in range(1, len(filtered_exon_pos[prot_id])):
        add_hang = fasta.get_seq(
            chrom,
            filtered_exon_pos[prot_id][n][0],
            filtered_exon_pos[prot_id][n][1],
        )
        left_hang = str(add_hang) + left_hang  # Prepend add_hang to left_hang
    # generate the overlapping windows
    left_hang = reverse_complement(left_hang)
    overlap = generate_windows(left_hang)
    return overlap


def flanking_lower_positions(string):
    first_upper_index = None
    last_upper_index = None

    for i, char in enumerate(string):
        if char.isupper():
            if first_upper_index is None:
                first_upper_index = i
            last_upper_index = i

    left_count = (
        sum(1 for i in range(first_upper_index) if string[i].islower())
        if first_upper_index is not None
        else 0
    )
    right_count = (
        sum(1 for i in range(last_upper_index + 1, len(string)) if string[i].islower())
        if last_upper_index is not None
        else 0
    )

    return left_count, right_count


def reverse_complement(dna_sequence):
    complement_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    reverse_comp_seq = "".join(
        complement_dict.get(base, base) for base in reversed(dna_sequence)
    )
    return reverse_comp_seq


def check_missing_codon(codon: str) -> bool:
    return codon != "-"


def perform_codon_check(
    var_class: str, codon_ref: str, codon_var: str
) -> (list[int], str, str):
    if var_class == "insertion":
        codon_check = codon_ref
    elif var_class == "deletion":
        codon_check = codon_var

    if check_missing_codon(codon_check):
        left_frames = [44]
    else:
        left_frames = [44, 45, 46]
        if var_class == "insertion":
            codon_ref = ""
        elif var_class == "deletion":
            codon_var = ""
    return left_frames, codon_ref, codon_var


def stop_check(aa: str) -> bool:
    return re.search("\\*", aa)


def frames_or_not(check: bool) -> list[int]:
    return [44, 45, 46] if check else [44]


def get_sequences_indel(
    mut_info,
    exon_info,
    fasta,
):
    (
        codon_ref,
        codon_var,
        chrom,
        mut,
        gene,
        strand,
        var_class,
    ) = (
        mut_info["codon_ref"],
        mut_info["codon_var"],
        mut_info["chromosome"],
        mut_info["mutation_location"],
        mut_info["gene_id"],
        mut_info["strand"],
        mut_info["variant_class"],
    )
    if gene not in exon_info:
        return None, None
    filtered_exon_pos = filter_exon_pos(exon_info[gene], mut, mut_info, strand)
    ref_seq_dict, var_seq_dict = {}, {}
    left_frames, codon_ref, codon_var = perform_codon_check(
        var_class, codon_ref, codon_var
    )
    for prot_id in filtered_exon_pos:
        ref_seq_dict[prot_id], var_seq_dict[prot_id] = defaultdict(list), defaultdict(
            list
        )
        for i, frame in enumerate(left_frames):
            if var_class == "insertion":
                codon_start_pos = (
                    mut - codon_var.find(next(filter(str.isupper, codon_var)))
                    if strand == "1"
                    else mut + codon_var.find(next(filter(str.isupper, codon_var)))
                )
            elif var_class == "deletion":
                codon_start_pos = (
                    mut - codon_ref.find(next(filter(str.isupper, codon_ref)))
                    if strand == "1"
                    else mut + codon_ref.find(next(filter(str.isupper, codon_ref)))
                )
            exon_start_diff = (
                codon_start_pos - filtered_exon_pos[prot_id][0][0]
                if strand == "1"
                else filtered_exon_pos[prot_id][0][1] - codon_start_pos
            )
            if exon_start_diff >= 44:
                start_pos = (
                    codon_start_pos - frame
                    if strand == "1"
                    else codon_start_pos + frame
                )  # this should generate a sequence of the left that is 45 bp long if +ve strand and to the right if -ve strand
            else:
                start_pos = (
                    filtered_exon_pos[prot_id][0][0]
                    if strand == "1"
                    else filtered_exon_pos[prot_id][0][1]
                )
                # check that the position from start_pos to codon_start_pos needs to be divisible by 3 and round down if not.
                if strand == "1":
                    if (codon_start_pos - 1 - start_pos) % 3 != 0:
                        start_pos = (
                            codon_start_pos
                            - ((codon_start_pos - start_pos) // 3) * 3
                            + 1
                        )
                else:
                    if (start_pos - codon_start_pos) % 3 != 0:
                        start_pos = (
                            codon_start_pos + ((start_pos - codon_start_pos) // 3) * 3
                        )
            if strand == "1":
                left_flank = fasta.get_seq(chrom, start_pos, codon_start_pos)
                # right side for ref
                end_pos_ref_pad = FINAL_SEQ_LEN - (len(left_flank) + len(codon_ref))
                end_pos_var_pad = FINAL_SEQ_LEN - (len(left_flank) + len(codon_var))
                if var_class == "insertion":
                    codon_end_pos = codon_start_pos + len(codon_ref) + 1
                elif var_class == "deletion":
                    codon_end_pos = codon_start_pos + len(codon_var) + 1
                end_pos_ref = codon_end_pos + end_pos_ref_pad - 1
                end_pos_var = codon_end_pos + end_pos_var_pad - 1
                exon_end_diff_ref = filtered_exon_pos[prot_id][0][1] - end_pos_ref
                exon_end_diff_var = filtered_exon_pos[prot_id][0][1] - end_pos_var
            else:
                right_flank = fasta.get_seq(chrom, codon_start_pos + 1, start_pos)
                # right side for ref
                end_pos_ref_pad = FINAL_SEQ_LEN - (len(right_flank) + len(codon_ref))
                end_pos_var_pad = FINAL_SEQ_LEN - (len(right_flank) + len(codon_var))
                if var_class == "insertion":
                    codon_end_pos = codon_start_pos - len(codon_ref)
                elif var_class == "deletion":
                    codon_end_pos = codon_start_pos - len(codon_var)
                end_pos_ref = codon_end_pos - end_pos_ref_pad
                end_pos_var = codon_end_pos - end_pos_var_pad
                exon_end_diff_ref = end_pos_ref - filtered_exon_pos[prot_id][0][0]
                exon_end_diff_var = end_pos_var - filtered_exon_pos[prot_id][0][0]
            if exon_end_diff_ref < 0:  # that means we exceed the end of the first exon
                if strand == "1":
                    if (
                        codon_end_pos > filtered_exon_pos[prot_id][0][1]
                    ):  # this means it's like to be an intronic mutation so it's irrelevant.
                        return None, None
                    right_flank_ref, right_flank_var = "", ""
                    if len(filtered_exon_pos[prot_id]) > 1:
                        right_flank_1 = fasta.get_seq(
                            chrom, codon_end_pos, filtered_exon_pos[prot_id][0][1]
                        )  # because this will be a fixed sequence
                        end_pos_ref = (
                            filtered_exon_pos[prot_id][1][0]
                            + abs(exon_end_diff_ref)
                            - 1
                        )
                        end_pos_var = (
                            filtered_exon_pos[prot_id][1][0]
                            + abs(exon_end_diff_var)
                            - 1
                        )
                        right_flank_2_ref = fasta.get_seq(
                            chrom, filtered_exon_pos[prot_id][1][0], end_pos_ref
                        )
                        right_flank_2_var = fasta.get_seq(
                            chrom, filtered_exon_pos[prot_id][1][0], end_pos_var
                        )
                        for part in [right_flank_1, right_flank_2_ref]:
                            right_flank_ref += str(part)
                        for part in [right_flank_1, right_flank_2_var]:
                            right_flank_var += str(part)
                    else:
                        return None, None  # because there's no more exons
                else:
                    if (
                        codon_end_pos < filtered_exon_pos[prot_id][0][0]
                    ):  # this means it's like to be an intronic mutation so it's irrelevant.
                        return None, None
                    left_flank_ref, left_flank_var = "", ""
                    if len(filtered_exon_pos[prot_id]) > 1:
                        left_flank_1 = fasta.get_seq(
                            chrom, codon_end_pos, filtered_exon_pos[prot_id][0][0]
                        )  # because this will be a fixed sequence
                        end_pos_ref = (
                            filtered_exon_pos[prot_id][1][1]
                            - abs(exon_end_diff_ref)
                            + 1
                        )
                        end_pos_var = (
                            filtered_exon_pos[prot_id][1][1]
                            - abs(exon_end_diff_var)
                            + 1
                        )
                        left_flank_2_ref = fasta.get_seq(
                            chrom, end_pos_ref + 1, filtered_exon_pos[prot_id][1][1]
                        )
                        left_flank_2_var = fasta.get_seq(
                            chrom, end_pos_var + 1, filtered_exon_pos[prot_id][1][1]
                        )
                        for part in [left_flank_1, left_flank_2_ref]:
                            left_flank_ref = str(part) + left_flank_ref
                        for part in [left_flank_1, left_flank_2_var]:
                            left_flank_var = str(part) + left_flank_var
                    else:
                        return None, None  # because there's no more exons
            else:
                if strand == "1":
                    right_flank_ref = fasta.get_seq(chrom, codon_end_pos, end_pos_ref)
                    right_flank_var = fasta.get_seq(chrom, codon_end_pos, end_pos_var)
                else:
                    left_flank_ref = fasta.get_seq(
                        chrom, end_pos_ref, codon_end_pos - 1
                    )
                    left_flank_var = fasta.get_seq(
                        chrom, end_pos_var, codon_end_pos - 1
                    )

            if strand == "1":
                ref_seq1 = str(left_flank) + codon_ref.upper() + str(right_flank_ref)
                var_seq1 = str(left_flank) + codon_var.upper() + str(right_flank_var)
            else:
                ref_seq1 = (
                    str(left_flank_ref)
                    + reverse_complement(codon_ref.upper())
                    + str(right_flank)
                )
                var_seq1 = (
                    str(left_flank_var)
                    + reverse_complement(codon_var.upper())
                    + str(right_flank)
                )
            if len(ref_seq1) != 93:
                return None, None
            if strand == "-1":
                ref_seq_dict[prot_id][i].append(reverse_complement(ref_seq1))
                var_seq_dict[prot_id][i].append(reverse_complement(var_seq1))
            else:
                ref_seq_dict[prot_id][i].append(ref_seq1)
                var_seq_dict[prot_id][i].append(var_seq1)
            if strand == "1":
                # get right hang
                overlap_ref = get_right_hang(
                    first_seq=ref_seq1,
                    chrom=chrom,
                    prot_id=prot_id,
                    first_seq_end_loc=end_pos_ref,
                    filtered_exon_pos=filtered_exon_pos,
                    fasta=fasta,
                )
                overlap_var = get_right_hang(
                    first_seq=var_seq1,
                    chrom=chrom,
                    prot_id=prot_id,
                    first_seq_end_loc=end_pos_var,
                    filtered_exon_pos=filtered_exon_pos,
                    fasta=fasta,
                )
            else:
                # get left hang. it will reverse complement at the end automatically.
                overlap_ref = get_left_hang(
                    first_seq=ref_seq1,
                    chrom=chrom,
                    prot_id=prot_id,
                    first_seq_start_loc=end_pos_ref,
                    filtered_exon_pos=filtered_exon_pos,
                    fasta=fasta,
                )
                overlap_var = get_left_hang(
                    first_seq=var_seq1,
                    chrom=chrom,
                    prot_id=prot_id,
                    first_seq_start_loc=end_pos_var,
                    filtered_exon_pos=filtered_exon_pos,
                    fasta=fasta,
                )
            ref_seq_dict[prot_id][i] += overlap_ref
            var_seq_dict[prot_id][i] += overlap_var

    return ref_seq_dict, var_seq_dict


def get_sequences_substitution(
    mut_info,
    exon_info,
    fasta,
):
    (
        codon_ref,
        codon_var,
        aa_ref,
        aa_var,
        chrom,
        mut,
        gene,
        strand,
    ) = (
        mut_info["codon_ref"],
        mut_info["codon_var"],
        mut_info["amino_acid_ref"],
        mut_info["amino_acid_var"],
        mut_info["chromosome"],
        mut_info["mutation_location"],
        mut_info["gene_id"],
        mut_info["strand"],
    )
    if gene not in exon_info:
        return None, None
    filtered_exon_pos = filter_exon_pos(exon_info[gene], mut, mut_info, strand)
    ref_seq_dict, var_seq_dict = {}, {}
    left_frames = frames_or_not(aa_var)
    for prot_id in filtered_exon_pos:
        ref_seq_dict[prot_id], var_seq_dict[prot_id] = defaultdict(list), defaultdict(
            list
        )
        for i, frame in enumerate(left_frames):
            codon_start_pos = (
                mut - codon_ref.find(next(filter(str.isupper, codon_ref)))
                if strand == "1"
                else mut + codon_ref.find(next(filter(str.isupper, codon_ref)))
            )  # always ref as there shouldn't be any "-" in codon_ref
            # also the codon_ref.find number works regardless of strand direction because the codon in the VEP output is already reverse complemented if it's negative strand
            exon_start_diff = (
                codon_start_pos - filtered_exon_pos[prot_id][0][0]
                if strand == "1"
                else filtered_exon_pos[prot_id][0][1] - codon_start_pos
            )
            if exon_start_diff >= 44:
                start_pos = (
                    codon_start_pos - frame
                    if strand == "1"
                    else codon_start_pos + frame
                )
                # this should generate a sequence of the left that is 45 bp long if +ve strand and to the right if -ve strand
            else:
                start_pos = (
                    filtered_exon_pos[prot_id][0][0]
                    if strand == "1"
                    else filtered_exon_pos[prot_id][0][1]
                )
                # check that the position from start_pos to codon_start_pos needs to be divisible by 3 and round down if not.
                if strand == "1":
                    if (codon_start_pos - 1 - start_pos) % 3 != 0:
                        start_pos = (
                            codon_start_pos
                            - ((codon_start_pos - start_pos) // 3) * 3
                            + 1
                        )
                else:
                    if (start_pos - codon_start_pos) % 3 != 0:
                        start_pos = (
                            codon_start_pos + ((start_pos - codon_start_pos) // 3) * 3
                        )
            if strand == "1":
                left_flank = fasta.get_seq(chrom, start_pos, codon_start_pos)
                end_pos_ref_pad = FINAL_SEQ_LEN - (len(left_flank) + len(codon_ref))
                end_pos_var_pad = FINAL_SEQ_LEN - (len(left_flank) + len(codon_var))
                codon_end_pos = codon_start_pos + len(codon_ref) + 1
                end_pos_ref = codon_end_pos + end_pos_ref_pad - 1
                end_pos_var = codon_end_pos + end_pos_var_pad - 1
                exon_end_diff_ref = filtered_exon_pos[prot_id][0][1] - end_pos_ref
                exon_end_diff_var = filtered_exon_pos[prot_id][0][1] - end_pos_var
            else:
                right_flank = fasta.get_seq(chrom, codon_start_pos + 1, start_pos)
                end_pos_ref_pad = FINAL_SEQ_LEN - (len(right_flank) + len(codon_ref))
                end_pos_var_pad = FINAL_SEQ_LEN - (len(right_flank) + len(codon_var))
                codon_end_pos = codon_start_pos - len(codon_ref)
                end_pos_ref = codon_end_pos - end_pos_ref_pad
                end_pos_var = codon_end_pos - end_pos_var_pad
                exon_end_diff_ref = end_pos_ref - filtered_exon_pos[prot_id][0][0]
                exon_end_diff_var = end_pos_var - filtered_exon_pos[prot_id][0][0]
            if exon_end_diff_ref < 0:  # that means we exceed the end of the first exon
                if strand == "1":
                    if (
                        codon_end_pos > filtered_exon_pos[prot_id][0][1]
                    ):  # this means it's like to be an intronic mutation so it's irrelevant.
                        return None, None
                    right_flank_ref, right_flank_var = "", ""
                    if len(filtered_exon_pos[prot_id]) > 1:
                        right_flank_1 = fasta.get_seq(
                            chrom, codon_end_pos, filtered_exon_pos[prot_id][0][1]
                        )  # because this will be a fixed sequence
                        end_pos_ref = (
                            filtered_exon_pos[prot_id][1][0]
                            + abs(exon_end_diff_ref)
                            - 1
                        )
                        end_pos_var = (
                            filtered_exon_pos[prot_id][1][0]
                            + abs(exon_end_diff_var)
                            - 1
                        )
                        right_flank_2_ref = fasta.get_seq(
                            chrom, filtered_exon_pos[prot_id][1][0], end_pos_ref
                        )
                        right_flank_2_var = fasta.get_seq(
                            chrom, filtered_exon_pos[prot_id][1][0], end_pos_var
                        )
                        for part in [right_flank_1, right_flank_2_ref]:
                            right_flank_ref += str(part)
                        for part in [right_flank_1, right_flank_2_var]:
                            right_flank_var += str(part)
                    else:
                        return None, None  # because there's no more exons
                else:
                    if (
                        codon_end_pos < filtered_exon_pos[prot_id][0][0]
                    ):  # this means it's like to be an intronic mutation so it's irrelevant.
                        return None, None
                    left_flank_ref, left_flank_var = "", ""
                    if len(filtered_exon_pos[prot_id]) > 1:
                        left_flank_1 = fasta.get_seq(
                            chrom, filtered_exon_pos[prot_id][0][0], codon_end_pos
                        )  # because this will be a fixed sequence
                        end_pos_ref = (
                            filtered_exon_pos[prot_id][1][1]
                            - abs(exon_end_diff_ref)
                            + 1
                        )
                        end_pos_var = (
                            filtered_exon_pos[prot_id][1][1]
                            - abs(exon_end_diff_var)
                            + 1
                        )
                        left_flank_2_ref = fasta.get_seq(
                            chrom, end_pos_ref + 1, filtered_exon_pos[prot_id][1][1]
                        )
                        left_flank_2_var = fasta.get_seq(
                            chrom, end_pos_var + 1, filtered_exon_pos[prot_id][1][1]
                        )
                        for part in [left_flank_1, left_flank_2_ref]:
                            left_flank_ref = str(part) + left_flank_ref
                        for part in [left_flank_1, left_flank_2_var]:
                            left_flank_var = str(part) + left_flank_var
                    else:
                        return None, None  # because there's no more exons
            else:
                if strand == "1":
                    right_flank_ref = fasta.get_seq(chrom, codon_end_pos, end_pos_ref)
                    right_flank_var = fasta.get_seq(chrom, codon_end_pos, end_pos_var)
                else:
                    left_flank_ref = fasta.get_seq(
                        chrom, end_pos_ref, codon_end_pos - 1
                    )
                    left_flank_var = fasta.get_seq(
                        chrom, end_pos_var, codon_end_pos - 1
                    )

            if strand == "1":
                ref_seq1 = str(left_flank) + codon_ref.upper() + str(right_flank_ref)
                var_seq1 = str(left_flank) + codon_var.upper() + str(right_flank_var)
            else:
                ref_seq1 = (
                    str(left_flank_ref)
                    + reverse_complement(codon_ref.upper())
                    + str(right_flank)
                )
                var_seq1 = (
                    str(left_flank_var)
                    + reverse_complement(codon_var.upper())
                    + str(right_flank)
                )
            if len(ref_seq1) != 93:
                return None, None
            if strand == "-1":
                ref_seq_dict[prot_id][i].append(reverse_complement(ref_seq1))
                var_seq_dict[prot_id][i].append(reverse_complement(var_seq1))
            else:
                ref_seq_dict[prot_id][i].append(ref_seq1)
                var_seq_dict[prot_id][i].append(var_seq1)
            if stop_check(aa_var):
                if strand == "1":
                    # get right hang
                    overlap_ref = get_right_hang(
                        first_seq=ref_seq1,
                        chrom=chrom,
                        prot_id=prot_id,
                        first_seq_end_loc=end_pos_ref,
                        filtered_exon_pos=filtered_exon_pos,
                        fasta=fasta,
                    )
                    overlap_var = get_right_hang(
                        first_seq=var_seq1,
                        chrom=chrom,
                        prot_id=prot_id,
                        first_seq_end_loc=end_pos_var,
                        filtered_exon_pos=filtered_exon_pos,
                        fasta=fasta,
                    )
                else:
                    # get left hang. it will reverse complement at the end automatically.
                    overlap_ref = get_left_hang(
                        first_seq=ref_seq1,
                        chrom=chrom,
                        prot_id=prot_id,
                        first_seq_start_loc=end_pos_ref,
                        filtered_exon_pos=filtered_exon_pos,
                        fasta=fasta,
                    )
                    overlap_var = get_left_hang(
                        first_seq=var_seq1,
                        chrom=chrom,
                        prot_id=prot_id,
                        first_seq_start_loc=end_pos_var,
                        filtered_exon_pos=filtered_exon_pos,
                        fasta=fasta,
                    )
                ref_seq_dict[prot_id][i] += overlap_ref
                var_seq_dict[prot_id][i] += overlap_var

    return ref_seq_dict, var_seq_dict
