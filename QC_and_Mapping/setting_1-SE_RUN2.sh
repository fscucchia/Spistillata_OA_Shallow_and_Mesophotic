in1=design/mapping_STAR.csv
out1=design/mapping_STAR.sh

python mapper_write-SE-or-PE-list_12.py --input_table $in1 --input_table_separator "," --output_script $out1 --SE
# --concatenate_partfiles --output_dir $outdir