### Concatenate files
- [`setting_1-SE.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/QC_and_Mapping/setting_1-SE.sh) that calls for [`mapper_write-SE-or_PE-list_12.py`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/QC_and_Mapping/mapper_write-SE-or-PE-list_12.py) and uses the design [`design/design1.csv`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/metadata/design1.csv).   
- [`mapper_write-SE-or_PE-list_12.py`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/QC_and_Mapping/mapper_write-SE-or-PE-list_12.py) creates the script templates [`design1.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/metadata/design1.sh) and [`design1.concatenated.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/metadata/design1.concatenated.sh).       
- [`design1.concatenated.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/metadata/design1.concatenated.sh) generates the concatenated files.          

### Quality filtering
- [`fastq-filter-SE_Conda_1_RUN1.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/QC_and_Mapping/fastq-filter-SE_Conda_1_RUN1.sh) uses the design file [`design1.concatenated.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/metadata/design1.concatenated.sh) to perform quality filtering (run starting from argument 0, then 1,2,3 and 5).        
- argument 2 calls for [`fastq-filter_job_5.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/QC_and_Mapping/fastq-filter_job_5.sh)(this script has all the cutadapt and trimmomatics commands) and removes the adapters using the file [`adapters4d.fa`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/QC_and_Mapping/adapters4d.fa).  

### Mapping
- [`setting_1-SE_RUN2.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/QC_and_Mapping/setting_1-SE_RUN2.sh) uses the design [`mapping_STAR.csv`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/metadata/mapping_STAR.csv) to create [`mapping_STAR.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/metadata/mapping_STAR.sh).
- [`star-conda_3_RUN1.sh`] uses the design [`mapping_STAR.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/metadata/mapping_STAR.sh) to map the filtered reads to _S. pistillata_ genome.


