function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }
################

# $1 is preserved for sbatch, while $2 and $3 are the steps
start=$2
end=$3
vcf1=$4
g1=$5
annotationDB=$6
cond1=$7
cond1mame=$8
#min_af=$9
#max_af=${10}
vcfRefDB=${9}
picard=${10} #/EV1/bioinfo/programs/picard/picard-tools-1.96
GATK=${11} #/EV1/bioinfo/programs/GATK3.1-1
java7=${12} #/EV1/bioinfo/programs/jdk1.7.0_60/bin/java
snpEff=${13} #/Home/programs/snpEff4.0E #snpEff=/EV1/bioinfo/programs/snpEff
gtf1=${14}

#################################################################################

if [ -e $vcf1 ]; then
	f1=$(echo $vcf1 | sed 's/.vcf//')
	if [ $start -le 0 ] && [ $end -ge 0 ]; then
		echo "vcf= "$vcf1
		echo "f1=  "$f1
		echo $annotationDB
		echo $cond1
		echo $cond1mame
	fi
	if [ $start -le 1 ] && [ $end -ge 1 ]; then
		cp $vcf1 $vcf1".backup"
		assert_
	fi
	if [ $start -le 2 ] && [ $end -ge 2 ]; then
		$java7 -jar $GATK"/"GenomeAnalysisTK.jar \
			-R $g1 \
			-T SelectVariants \
			--variant "$vcf1" \
			-o "$f1.indels.vcf" \
			-selectType INDEL
		assert_
		$java7 -jar $GATK"/"GenomeAnalysisTK.jar \
			-R $g1 \
			-T SelectVariants \
			--variant "$vcf1" \
			-o "$f1.snps.vcf" \
			-selectType SNP
		assert_
	fi
	if [ $start -le 3 ] && [ $end -ge 3 ]; then
		# filter = remove those variants which satisfy the condition
		echo $cond1 > "$f1.$cond1mame.txt"
		assert_
		$java7 -jar $GATK"/"GenomeAnalysisTK.jar \
			-T VariantFiltration \
			-R $g1 \
			-V "$f1.indels.vcf" \
			-window 35 \
			-cluster 3 \
			-filterName FS \
			-filter "FS > 30.0" \
			-filterName QD \
			-filter "QD < 2.0" \
			-filterName "$cond1mame" \
			-filter "$cond1" \
			-o "$f1.$cond1mame.indels.vcf"
		assert_
		grep 'PASS\|^#' "$f1.$cond1mame.indels.vcf" > "$f1.$cond1mame.indels.pass.vcf"
        assert_
        #perl -ne 'use List::Util qw[min max]; $x=$_; chomp($x); $x =~ m/AF=([0-9\.\,]+);/; @af=split(/,/,$1); if(((max(@af)>'$min_af') and (max(@af)<'$max_af'))or($x =~ m/^#.*/)){print $_}'
		#grep 'PASS\|^#' "$f1.$cond1mame.indels.vcf" > "$f1.$cond1mame.indels.badAF.vcf"
        #assert_
        # perl -ne 'use List::Util qw[min max]; $x=$_; chomp($x); $x =~ m/AF=([0-9\.\,]+);/; @af=split(/,/,$1); if(((max(@af)>'$min_af') and (max(@af)<'$max_af'))or($x =~ m/^#.*/)){}else{print $_}'
		$java7 -jar $GATK"/"GenomeAnalysisTK.jar \
			-T VariantFiltration \
			-R $g1 \
			-V "$f1.snps.vcf" \
			-window 35 \
			-cluster 3 \
			-filterName FS \
			-filter "FS > 30.0" \
			-filterName QD \
			-filter "QD < 2.0" \
			-filterName "$cond1mame" \
			-filter "$cond1" \
			-o "$f1.$cond1mame.snps.vcf"
		assert_
		grep 'PASS\|^#' "$f1.$cond1mame.snps.vcf" > "$f1.$cond1mame.snps.pass.vcf"
        assert_
        # perl -ne 'use List::Util qw[min max]; $x=$_; chomp($x); $x =~ m/AF=([0-9\.\,]+);/; @af=split(/,/,$1); if(((max(@af)>'$min_af') and (max(@af)<'$max_af'))or($x =~ m/^#.*/)){print $_}'
		#grep 'PASS\|^#' "$f1.$cond1mame.snps.vcf" > "$f1.$cond1mame.snps.badAF.vcf"
        #assert_
        # perl -ne 'use List::Util qw[min max]; $x=$_; chomp($x); $x =~ m/AF=([0-9\.\,]+);/; @af=split(/,/,$1); if(((max(@af)>'$min_af') and (max(@af)<'$max_af'))or($x =~ m/^#.*/)){}else{print $_}'
	fi
	if [ $start -le 4 ] && [ $end -ge 4 ]; then
        #. "$CONDA"
		#conda activate snpeff_setting1
		
        $java7 -Xmx128g -jar $snpEff/snpEff/snpEff.jar \
            "$annotationDB" \
			-c $snpEff/snpEff/snpEff.config \
            -v \
			-s "$f1.$cond1mame.indels.pass.snpEff.html" \
			"$f1.$cond1mame.indels.pass.vcf" > \
			"$f1.$cond1mame.indels.pass.snpEff.vcf"
		assert_
        
		$java7 -Xmx128g -jar $snpEff/snpEff/snpEff.jar \
            "$annotationDB" \
			-c $snpEff/snpEff/snpEff.config \
            -v \
			-s "$f1.$cond1mame.snps.pass.snpEff.html" \
			"$f1.$cond1mame.snps.pass.vcf" > \
			"$f1.$cond1mame.snps.pass.snpEff.vcf"
		assert_
	fi
	#if [ $start -le 5 ] && [ $end -ge 5 ]; then
	#	bedtools intersect -wb -a "$f1.$cond1mame.snps.pass.vcf"   -b $gtf1 > "$f1.$cond1mame.snps.pass.annotate.bed"
	#	assert_
	#	bedtools intersect -wb -a "$f1.$cond1mame.indels.pass.vcf" -b $gtf1 > "$f1.$cond1mame.indels.pass.annotate.bed"
	#	assert_
	#	cat "$f1.$cond1mame.snps.pass.annotate.bed"   | grep 'CDS'  |  sort -t $'\t' -k1,1 -k2,2n | awk 'BEGIN {scaffold=""; pos=""} scaffold != $1 || pos != $2 {print; scaffold=$1; pos=$2}' > "$f1.$cond1mame.snps.pass.annotate.CDS.bed"
	#	assert_
	#	cat "$f1.$cond1mame.indels.pass.annotate.bed" | grep 'CDS'  |  sort -t $'\t' -k1,1 -k2,2n | awk 'BEGIN {scaffold=""; pos=""} scaffold != $1 || pos != $2 {print; scaffold=$1; pos=$2}' > "$f1.$cond1mame.indels.pass.annotate.CDS.bed"
	#	assert_
	#	cat "$f1.$cond1mame.snps.pass.annotate.bed"   | grep 'exon' |  sort -t $'\t' -k1,1 -k2,2n | awk 'BEGIN {scaffold=""; pos=""} scaffold != $1 || pos != $2 {print; scaffold=$1; pos=$2}' > "$f1.$cond1mame.snps.pass.annotate.exon.bed"
	#	assert_
	#	cat "$f1.$cond1mame.indels.pass.annotate.bed" | grep 'exon' |  sort -t $'\t' -k1,1 -k2,2n | awk 'BEGIN {scaffold=""; pos=""} scaffold != $1 || pos != $2 {print; scaffold=$1; pos=$2}' > "$f1.$cond1mame.indels.pass.annotate.exon.bed"
	#	assert_
	#fi
	#if [ $start -le 6 ] && [ $end -ge 6 ]; then
	#	# snps:
	#	# annotation using existing snps in a vcf file
	#	$java7 -jar $snpEff/SnpSift.jar annotate "$vcfRefDB" "$f1.$cond1mame.snps.pass.snpEff.vcf" > "$f1.$cond1mame.snps.pass.snpSift.vcf"
	#	assert_
	#	# vcf to table including snpEff data
	#	$java7 -Xmx128g -jar $snpEff/SnpSift.jar extractFields -s "|" "$f1.$cond1mame.snps.pass.snpSift.vcf" \
	#		"CHROM" "POS" "REF" "ALT" "QUAL" "AF" "DP" "MQ" \
	#		"EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].TRID" \
	#		"GEN[*].GT" "GEN[*].GT" "GEN[*].AD" "GEN[*].DP" \
	#		> "$f1.$cond1mame.snps.pass.snpSift.txt"
	#	assert_
	#	#$java7 -jar $snpEff/SnpSift.jar tstv hom "$f1.$cond1mame.snps.pass.snpEff.vcf" > "$f1.$cond1mame.snps.pass.snpSift.tstv.txt"
	#	#assert_
	#	
	#	# indels:
	#	# annotation using existing snps in a vcf file
	#	$java7 -jar $snpEff/SnpSift.jar annotate "$vcfRefDB" "$f1.$cond1mame.indels.pass.snpEff.vcf" > "$f1.$cond1mame.indels.pass.snpSift.vcf"
	#	assert_
	#	# vcf to table including snpEff data
	#	$java7 -Xmx128g -jar $snpEff/SnpSift.jar extractFields -s '|' "$f1.$cond1mame.indels.pass.snpSift.vcf" \
	#		"CHROM" "POS" "REF" "ALT" "QUAL" "AF" "DP" "MQ" \
	#		"EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].TRID" \
	#		"GEN[*].GT" "GEN[*].GT" "GEN[*].AD" "GEN[*].DP" \
	#		> "$f1.$cond1mame.indels.pass.snpSift.txt"
	#	assert_
	#fi
	if [ $start -le 7 ] && [ $end -ge 7 ]; then
		# http://www.broadinstitute.org/gatk/guide/article?id=1268
		# QUAL = -10*log10(1-p) , where p is the probability of true call (10 indicates a 1/10 chance of error, 100 indicates 1/00 ...)
		# -raw = default is false
		
		$java7 -jar $GATK"/"GenomeAnalysisTK.jar \
		     -R $g1 \
		     -T VariantsToTable \
		     -V "$f1.$cond1mame.snps.pass.vcf" \
		     -F CHROM -F POS -F FILTER -F REF -F ALT -F AC -F AN -F AF -F QUAL -F DP -F MQ -F TRANSITION \
		     -o "$f1.$cond1mame.snps.pass.txt"
		assert_
		
		$java7 -jar $GATK"/"GenomeAnalysisTK.jar \
		     -R $g1 \
		     -T VariantsToTable \
		     -V "$f1.$cond1mame.indels.pass.vcf" \
		     -F CHROM -F POS -F FILTER -F REF -F ALT -F AC -F AN -F AF -F QUAL -F DP -F MQ -F TRANSITION  \
		     -o "$f1.$cond1mame.indels.pass.txt"
		assert_
	fi
    if [ $start -le 8 ] && [ $end -ge 8 ]; then
        cat "$f1.$cond1mame.snps.pass.txt" | perl -ne 'chomp; @x=split("\t",$_); printf("%s_%s\t%s\n",$x[0],$x[1],$_)' > "$f1.$cond1mame.snps.pass.2.txt"
        assert_
        cat "$f1.$cond1mame.indels.pass.txt" | perl -ne 'chomp; @x=split("\t",$_); printf("%s_%s\t%s\n",$x[0],$x[1],$_)' > "$f1.$cond1mame.indels.pass.2.txt"
        assert_
    fi
fi


