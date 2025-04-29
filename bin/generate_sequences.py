#!/usr/bin/env python
import argparse

from pathlib import Path
from pyfaidx import Fasta


from minigene_design._utils import (
    MINIGENE_SEQ_LEN,
    KEEP_FROM_LINES,
    read_and_filter,
    extract_exon_info,
    create_minigenes,
    prepare_output,
)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-r",
        "--ref",
        type=str,
        help="path to the reference genome in FASTA format.",
    )
    parser.add_argument(
        "-g",
        "--gtf",
        type=str,
        help="path to the GTF file with exon information.",
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help="path to the input vep table.",
    )
    parser.add_argument(
        "-o",
        "--out_dir",
        type=str,
        help="output directory.",
    )
    parser.add_argument(
        "-s",
        "--sample_id",
        type=str,
        default="sample",
        help="sample id.",
    )
    parser.add_argument(
        "-k",
        "--keep_from_lines",
        type=int,
        default=KEEP_FROM_LINES,
        help="number of lines to keep from the input vep table.",
    )
    parser.add_argument(
        "-l",
        "--return_length",
        type=int,
        default=MINIGENE_SEQ_LEN,
        help="length of the minigene sequence to return.",
    )
    parser.add_argument(
        "-c",
        "--correct_bbsl",
        action="store_true",
        help="correct the BBSL for the minigene.",
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    refgen = Fasta(args.ref)
    exon_info = extract_exon_info(args.gtf)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = read_and_filter(args.input, keep_from_lines=args.keep_from_lines)
    result = {}
    result = create_minigenes(df=df, sample=args.sample_id, exon_info=exon_info, fasta=refgen, out_dict=result)
    ref_df, var_df = prepare_output(result, correct_bbsl=args.correct_bbsl)
    ref_df.to_csv(out_dir / f"{args.sample_id}_ref.csv")
    var_df.to_csv(out_dir / f"{args.sample_id}_var.csv")


if __name__ == "__main__":
    main()
