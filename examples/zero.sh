# install
cd /path/to/minigene_design
pip install .

# run
python bin/generate_sequences.py \
    -r /path/to/resources/GRCh38.d1.vd1.fa \
    -g /path/to/resources/gencode.v33.primary_assembly.annotation.gtf \
    -i test/test_hc_vep.txt \
    -o test_out \
    --keep_from_lines 93 \
    --return_length 81
