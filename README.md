# Genetic and physiological traits conferring tolerance to ocean acidification in mesophotic corals  [![DOI:10.1111/gcb.15812](http://img.shields.io/badge/DOI-10.1111/gcb.15812-33C4FF.svg)](https://doi.org/10.1111/gcb.15812)
Federica Scucchia<sup>1,2</sup>, Assaf Malik<sup>1</sup>, Hollie M. Putnam<sup>3</sup>, Tali Mass<sup>1,4</sup>

<sup>1</sup> Department of Marine Biology, Leon H. Charney school of Marine Sciences, University of Haifa, Haifa, 3498838, Israel                                                                               
<sup>2</sup> The Interuniversity Institute of Marine Sciences, Eilat 88103, Israel                                                                                                             
<sup>3</sup> Department of Biological Sciences, University of Rhode Island, Kingston, RI02881, United States of America                                                                      
<sup>4</sup> Morris Kahn Marine Research Station, The Leon H. Charney School of Marine Sciences, University of Haifa, Sdot Yam, Israel

![pic](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/blob/main/media/Shallow_and_Mesophotic_pictures.jpg?raw=true)

This electronic notebook provides the scripts employed to analyze _Stylophora pistillata_ gene expression dynamics across shallow (5 m) and mesophotic (45 m) reefs under ambient (8.2 pH), intermediate-low (7.8 pH), and low (7.6 pH) pH conditions. This analysis characterizes the transcriptomic response of adult _S. pistillata_ colonies from the Red Sea to predicted future ocean acidification conditions.

### RNA-Seq reads quality filtering and mapping

Total RNA was extracted from shallow and mesophotic _S. pistillata_ fragments collected in Eilat (Red Sea). RNA-seq libraries were prepared using an in-house protocol at the Weizmann Institute of Science (Israel). 

**[Raw sequence data - (NCBI SRA)](https://dataview.ncbi.nlm.nih.gov/object/PRJNA701170?reviewer=2sdh5ejluhr6na11607otr0c7i)** - Raw Illumina sequence data (fastq format) for the 6 coral colonies succesfully sequenced in this study (accession code [PRJNA701170](https://dataview.ncbi.nlm.nih.gov/object/PRJNA701170?reviewer=2sdh5ejluhr6na11607otr0c7i)).

**[Quality filtering and mapping](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/tree/main/QC_and_Mapping)** - Details the processing and analyses of the _S. pistillata_ transcriptome sequencing data. RNA-Seq reads processing included adapter trimming using Cutadapt v1.15 ([Martin, 2011](https://doi.org/10.14806/ej.17.1.200)) and quality filtering using Trimmomatic v0.3 ([Bolger et al., 2014](https://doi.org/10.1093/bioinformatics/btu170)). Reads were aligned to the host genome assembly [NCBI GCA_002571385.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_002571385.1/) using STAR v2.5 ([Dobin et al., 2013](https://doi.org/10.1093/bioinformatics/bts635); [Voolstra et al., 2017](https://doi.org/10.1038/s41598-017-17484-x)). 

### Species identification

**[Coral host and symbiont species identification](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/tree/main/Species_identification)** - High quality reads were blasted using Diamond v2.0.7 ([Buchfink et al., 2015](https://doi.org/10.1038/nmeth.3176)) against the [NCBI](https://www.ncbi.nlm.nih.gov/) and [Reefgenomics](http://reefgenomics.org/) genome-based proteomes databases of Symbiodiniaceae species _Symbiodinium microadraticum_, _Cladocopium goreaui_, and _Fugacium kawagutii_ (formerly _Symbiodinium_ spp. clades A, C1, and F, respectively ([LaJeunesse et al., 2018](https://doi.org/10.1016/j.cub.2018.07.008))), as well as to the databases of Cnidaria, selected stramenopiles/alveolates/Rhizaria and Metazoa.

### Differential expression

**[Coral host differential expression](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/tree/main/DE)** - DE analysis was conducted using Bioconductor DEseq2 v1.30.1 ([Love et al., 2014](https://doi.org/10.1186/s13059-014-0550-8)) by a) analyzing mesophotic and shallow samples separately, considering a single factor (pH) and three levels (8.2, 7.8, 7.6), and b) analyzing samples from each pH separately, considering a single factor (depth) and two levels (shallow, mesophotic).

### Functional enrichment

**[Coral host functional enrichment](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/tree/main/Enrichment)** - Genes biological terms were assigned based on _S. pistillata_ [Uniprot](https://www.uniprot.org/), [KEGG](https://www.kegg.jp/) and Trinotate annotations ([Bryant et al., 2017](https://doi.org/10.1016/j.celrep.2016.12.063)). Enrichment analysis was conducted in Bioconductor GOSeq v1.42.0, as described previously ([Malik et al., 2020](https://doi.org/10.1016/j.actbio.2020.01.010)). 

### SNPs characterization

**[Coral host SNPs characterization](https://github.com/fscucchia/Spistillata_OA_Shallow_and_Mesophotic/tree/main/SNPs)** - SNPs analysis was conducted using the recommended RNA-Seq SNPs practice of the [Broad Institute](https://gatk.broadinstitute.org). Alignment biases were removed from the STAR-aligned reads, variant calling was performed using the GATK HaplotypeCaller v3.5 ([McKenna et al., 2010](https://doi.org/10.1101/gr.107524.110)), and finally variant filtration and annotation was carried out following ([DePristo et al., 2011](https://doi.org/10.1038/ng.806)). 
Genetic differentiation among depths was assessed using Identity By State (IBS) analysis, which was conducted using [SNPRelate package v1.24.0](https://github.com/zhengxwen/SNPRelate) in R. The degree of genetic differentiation was also determined by estimating the fixation index (_Fst_)([Weir & Cockerham, 1984](https://doi.org/10.1111/j.1558-5646.1984.tb05657.x)) using the [HIERFSTAT package v0.5-7](https://cran.r-project.org/web/packages/hierfstat/index.html) in R, between samples located at the same depth and between samples at different depths. 
