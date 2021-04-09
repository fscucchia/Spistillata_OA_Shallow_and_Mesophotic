in1=design/design1.csv
out1=design/design1.sh
outdir="/Home/concat"

python mapper_write-SE-or-PE-list_12.py --input_table $in1 --input_table_separator "," --output_script $out1 --concatenate_partfiles --output_dir $outdir --SE
