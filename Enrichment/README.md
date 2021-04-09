### Functional enrichment 
- `goseq-foreground-background_4_RUN1.R` (GOSeq tables settings), take the result os DESeq and makes tables (this script is for the comparisons mesophotic vs shallow within the same pH)         x
- `goseq-foreground-background_4_RUN2.R` (GOSeq tables settings), take the result os DESeq and makes tables (this script is for the comparisons among pH conditions within the same depth)       x
- run using the command *bash* `goseq_enrichment_RUN1.sh` (GOSeq stylophora enrichment based on stylophora KEGG , GO-Uniprot, GO-Trinotate databases) arguments 0 to 2, that uses the DE_file='/enrichement/shallow.comparison.txt' (created by goseq-foreground-background_4_RUN1/2.R) and that calls for `goseq_enrichment.R` (GOSeq call)
