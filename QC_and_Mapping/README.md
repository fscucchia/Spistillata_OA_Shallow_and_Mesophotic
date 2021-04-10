### Concatenate files
- [`setting_1-SE.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/QC_and_Mapping/setting_1-SE.sh) that calls for [`mapper_write-SE-or_PE-list_12.py`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/QC_and_Mapping/mapper_write-SE-or-PE-list_12.py) and uses the design [`design/design1.csv`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/metadata/design1.csv).   
- [`mapper_write-SE-or_PE-list_12.py`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/QC_and_Mapping/mapper_write-SE-or-PE-list_12.py) creates the script templates [`design1.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/metadata/design1.sh) and [`design1.concatenated.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/metadata/design1.concatenated.sh).       
- [`design1.concatenated.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/metadata/design1.concatenated.sh) generates the concatenated files.          

### Quality filtering
- `fastq-filter-SE_Conda_1_RUN1.sh` that uses the design file `design1.concatenated.sh`          x
- run using the command *bash* `fastq-filter-SE_Conda_1_RUN1.sh` starting from argument 0, then 1,2,3 and 5 (step 5 is multiqc)            x
- argument 2 of `fastq-filter-SE_Conda_1_RUN1.sh` calls for `fastq-filter_job_5.sh`(this script has all the cutadapt and trimmomatics commands) and removes the adapters using the file `adapters4d.fa` and performs quality filtering             x

### Mapping
- `setting_1-SE_RUN2.sh` uses the design `mapping_STAR.csv` and creates `mapping_STAR.sh`      x
- run using the command *bash* `setting_1-SE_RUN2.sh`(it will call for `mapper_write-SE-or_PE-list_12.py`)          x
- `star-conda_3_RUN1.sh` (uses the design `mapping_STAR.sh`, and Stylophora genome) starting from argument 1, then 2,3 and 5 (step 5 is multiqc)      x
