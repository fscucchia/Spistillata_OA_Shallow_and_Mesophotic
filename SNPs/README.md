### Processing of STAR-aligned reads 
- [`gatk21_threads_RUN1.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/SNPs/gatk21_threads_RUN1.sh) calls for [`gatpk21_pipeline.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/SNPs/gatk21_pipeline.sh) and [`gatpk21_vcf-processing.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/SNPs/gatk21_vcf-processing.sh) to remove alignment biases from the STAR-aligned reads, and to perform variant calling (using the GATK HaplotypeCaller), variant filtration and annotation.

### Identity By State (IBS) analysis
- [`joint-snps4_RUN1.R`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/SNPs/joint-snps4_RUN1.R) performs the IBS analysis using the script [`bams-mos-depth_1.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/SNPs/bams-mos-depth_1.sh) and the design [`design2.csv`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/metadata/design2.csv).

###  Fst estimates
- [`Fst_statistics.R`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/SNPs/Fst_statistics.R) takes the output of [`joint-snps4_RUN1.R`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/SNPs/joint-snps4_RUN1.R) to prepare the input table for hierfstat (a [data frame](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/SNPs/SNPs_locus_table.txt) with genotype of the individuals -one locus per column) and calculates Fst statistics.



