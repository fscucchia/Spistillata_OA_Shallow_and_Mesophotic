in1=design/design1.csv #https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/metadata/design1.csv
out1=design/design1.sh #https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/metadata/design1.sh
outdir="/Home/concat"

# mapper_write-SE-or-PE-list_12.py (https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/QC_and_Mapping/mapper_write-SE-or-PE-list_12.py)
python mapper_write-SE-or-PE-list_12.py --input_table $in1 --input_table_separator "," --output_script $out1 --concatenate_partfiles --output_dir $outdir --SE
