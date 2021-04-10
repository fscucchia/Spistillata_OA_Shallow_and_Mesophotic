### Processing of STAR-aligned reads 
- [`gatk21_threads_RUN1.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/SNPs/gatk21_threads_RUN1.sh) calls for [`gatpk21_pipeline.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/SNPs/gatk21_pipeline.sh) and [`gatpk21_vcf-processing.sh`](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/SNPs/gatk21_vcf-processing.sh) to remove alignment biases from the STAR-aligned reads, and to perform variant calling (using the GATK HaplotypeCaller), variant filtration and annotation.

### Identity By State (IBS) analysis:

- joint-snps4_RUN1.R
- bams-mos-depth_1.sh
- design3.csv #metadata

###  Fst estimates:

- Fst_statistics.R
- SNPs_locus_table #a data frame with the genotype of the individuals -one locus per column

