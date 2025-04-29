# minigene_design

## Installation

```bash
git clone https://github.com/tuonglab/minigene_design.git
cd minigene_design
pip install .
```

## Quick start

```bash
python bin/generate_sequences.py \
    -r /path/to/resources/GRCh38.d1.vd1.fa \
    -g /path/to/resources/gencode.v33.primary_assembly.annotation.gtf \
    -i test/test_hc_vep.txt \
    -o test_out
```

```bash
python bin/generate_sequences.py -h
usage: generate_sequences.py [-h] [-r REF] [-g GTF] [-i INPUT] [-o OUT_DIR] [-s SAMPLE_ID] [-k KEEP_FROM_LINES] [-l RETURN_LENGTH] [-c]

options:
  -h, --help            show this help message and exit
  -r, --ref REF         path to the reference genome in FASTA format.
  -g, --gtf GTF         path to the GTF file with exon information.
  -i, --input INPUT     path to the input vep table.
  -o, --out_dir OUT_DIR
                        output directory.
  -s, --sample_id SAMPLE_ID
                        sample id.
  -k, --keep_from_lines KEEP_FROM_LINES
                        number of lines to keep from the input vep table.
  -l, --return_length RETURN_LENGTH
                        length of the minigene sequence to return.
  -c, --correct_bbsl    correct the BBSL for the minigene.
```
