
#################################################################################

preserved=$1
start=$2
end=$3
f1=$4
g1=$5
outdir1=$6
whichGenotyper=$7
name1=$8
phred=$9
bamVersion=${10} # bamVersion = 5 (GATK original setting), or 6 (after removing all multireads with bamutils) 
ploidy=${11} 
vcfRefDB1=${12}
picard=${13} 
GATK=${14} 
java7=${15} 
snpEff=${16} 
vcfPrevSycle1=${17}
cycle=${18}

#################################################################################

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

echo "step="$start"-"$end

if [ $start -ge -1 ] && [ $end -le 12 ]; then
	echo "f1="$f1
	echo "name="$name1
	echo "outdir="$outdir1
	if [ $start -le 0 ] && [ $end -ge 0 ]; then # RUN THIS STEP ( step 0 ) ONLY ONCE !!!! 
		echo "test: ref = "$g1
		echo "test: out = "$outdir1
		
		echo "creating faidx for $g1"
		samtools faidx $g1
		assert_
		dict=$(echo $g1 | sed 's/.fa$\|.fasta$//')".dict"
		echo "dict="$dict
		echo "genome="$g1
		echo "creating dict ..."
		java -jar $picard CreateSequenceDictionary REFERENCE=$g1 OUTPUT=$dict
		assert_
	fi
	if [ $start -le 1 ] && [ $end -ge 1 ]; then # index bams
		samtools index $f1
		assert_
	fi
	if [ $start -le 2 ] && [ $end -ge 2 ]; then # make output dirs
		mkdir $outdir1"/"$name1
		assert_
		mkdir $outdir1"/"$name1"/TEMP"
		assert_
	fi
	if [ $start -le 3 ] && [ $end -ge 3 ]; then
		in1=$f1
		lib1=$name1
		java -jar $picard AddOrReplaceReadGroups \
			INPUT=$in1 \
			OUTPUT="$outdir1/$name1/$name1.2.bam" \
			SORT_ORDER=coordinate \
			RGID=$name1 \
			RGLB=$lib1 \
			RGPL='ILLUMINA' \
			RGPU='machine1' \
			RGSM=$name1  \
			TMP_DIR=$outdir1"/"$name1"/TEMP"
			#RGID='group1' \
		assert_
	fi
	if [ $start -le 4 ] && [ $end -ge 4 ]; then
		java -jar $picard MarkDuplicates \
			INPUT="$outdir1/$name1/$name1.2.bam" \
			OUTPUT="$outdir1/$name1/$name1.3.bam"  \
			CREATE_INDEX=true \
			VALIDATION_STRINGENCY=SILENT \
			METRICS_FILE="$outdir1/$name1/$name1.3.output.metrics" \
			TMP_DIR=$outdir1"/"$name1"/TEMP"
		assert_
	fi
	if [ $start -le 5 ] && [ $end -ge 5 ]; then
		if [ $phred -eq 33 ]; then
			$java7 -jar $GATK"/"GenomeAnalysisTK.jar \
				-T SplitNCigarReads \
				-R $g1 \
				-I "$outdir1/$name1/$name1.3.bam" \
				-o "$outdir1/$name1/$name1.4.bam" \
				-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 \
				-U ALLOW_N_CIGAR_READS
			assert_
		elif [ $phred -eq 64 ]; then
			$java7 -jar $GATK"/"GenomeAnalysisTK.jar \
				-T SplitNCigarReads \
				-R $g1 \
				-I "$outdir1/$name1/$name1.3.bam" \
				-o "$outdir1/$name1/$name1.4.bam" \
				-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 \
				-U ALLOW_N_CIGAR_READS \
				--fix_misencoded_quality_scores
			assert_
		else
			exit "unrecognized fastq quality format "$phred
		fi
	fi
	if [ $start -le 6 ] && [ $end -ge 6 ]; then
		$java7 -jar $GATK"/"GenomeAnalysisTK.jar \
    		-T RealignerTargetCreator \
			-R $g1 \
			-I "$outdir1/$name1/$name1.4.bam" \
			-o "$outdir1/$name1/$name1.4.target_intervals.list"
		assert_
		#-known gold_indels.vcf
	fi
	if [ $start -le 7 ] && [ $end -ge 7 ]; then
		$java7 -jar $GATK"/"GenomeAnalysisTK.jar \
		    -T IndelRealigner \
		    -R $g1 \
		    -I "$outdir1/$name1/$name1.4.bam" \
		    -targetIntervals "$outdir1/$name1/$name1.4.target_intervals.list" \
		    -o "$outdir1/$name1/$name1.5.bam"
		assert_
		#-known gold_indels.vcf
	fi	
	if [ $start -le 80 ] && [ $end -ge 80 ]; then
		echo "bam: ""$outdir1/$name1/$name1.5.bam"
		bamutils filter "$outdir1/$name1/$name1.5.bam" "$outdir1/$name1/$name1.6.bam" -failed "$outdir1/$name1/$name1.6.multimappers" -lt NH:i 2
		assert_
		samtools index "$outdir1/$name1/$name1.6.bam"
		assert_
	fi
	if [ $start -le 90 ] && [ $end -ge 90 ]; then
		b1=$bamVersion # 5 or 6
		# RECALIBRATION	USING VCF CREATED IN PREVIOUS CYCLE AND AN EXITING REFERENCE
		$java7 -jar $GATK"/"GenomeAnalysisTK.jar \
			   -T BaseRecalibrator \
			   -I "$outdir1/$name1/$name1.$b1.bam" \
			   -R $g1 \
		           -knownSites $vcfRefDB1 \
		           -knownSites $vcfPrevSycle1 \
			   -o "$outdir1/$name1/$name1.7.cycle$cycle.grp"
		assert_
		# PRINTING RECALIBRATED BAM FILE
		$java7 -jar $GATK"/"GenomeAnalysisTK.jar \
   			-T PrintReads \
   			-R $g1 \
   			-I "$outdir1/$name1/$name1.$b1.bam" \
   			-BQSR "$outdir1/$name1/$name1.7.cycle$cycle.grp" \
   			-o "$outdir1/$name1/$name1.7.bam"
		assert_
	fi
	if [ $start -le 10 ] && [ $end -ge 10 ]; then
		b1=$bamVersion # 5 or 6
		if [ $whichGenotyper -eq 0 ]; then # test
			echo "test name="$name1
			echo "test R="$g1
			if [ -f "$g1" ]; then echo 'OK'; else echo 'file not found'; fi
			echo "test I=""$outdir1/$name1/$name1.$b1.bam"
			if [ -f "$outdir1/$name1/$name1.$b1.bam" ]; then echo 'OK'; else echo 'file not found'; fi
			echo "test o=""$outdir1/$name1/$name1.$b1.vcf"
			echo "ploidy="$ploidy
			echo '------------------------------------'
		fi
		if [ $whichGenotyper -gt 0 ]; then
		
			echo 'HaplotypeCaller single sample discovery'
			# single sample discovery
			# -recoverDanglingHeads is deprecated (became default since GATK3.3, therefore removed)
			$java7 -jar $GATK"/"GenomeAnalysisTK.jar \
				-T HaplotypeCaller \
				-R $g1 \
				-I "$outdir1/$name1/$name1.$b1.bam" \
				-dontUseSoftClippedBases \
				-stand_call_conf 20.0 \
				-stand_emit_conf 20.0 \
				-o "$outdir1/$name1/$name1.$b1.vcf" \
				--sample_ploidy $ploidy
			assert_
			
			echo 'HaplotypeCaller GVCF for joint discovery'
			# for joint discovery (check if approved for for rna-seq)
			# -recoverDanglingHeads is deprecated (became default since GATK3.3, therefore removed)
			$java7 -jar $GATK"/"GenomeAnalysisTK.jar \
				-T HaplotypeCaller \
				-R $g1 \
				-I "$outdir1/$name1/$name1.$b1.bam" \
				-dontUseSoftClippedBases \
				-stand_call_conf 20.0 \
				-stand_emit_conf 20.0 \
				-o "$outdir1/$name1/$name1.$b1.g.vcf" \
				--emitRefConfidence GVCF \
				--sample_ploidy $ploidy
				# --variant_index_type LINEAR --variant_index_parameter 128000 \ # versions older than 3.4 require that
				# -o "$outdir1/$name1/$name1.$b1.gGVCF.vcf" \
			assert_
		fi
	fi
	if [ $start -le 11 ] && [ $end -ge 11 ]; then
		b1=$bamVersion # 5 or 6
		mkdir "$outdir1/$name1/filtered1"
		$java7 	-jar $GATK"/"GenomeAnalysisTK.jar \
			-T VariantFiltration \
			-R $g1 \
			-V "$outdir1/$name1/$name1.$b1.vcf" \
			-window 35 -cluster 3 \
			-filterName FS -filter "FS > 30.0" \
			-filterName QD -filter "QD < 2.0" \
			-o "$outdir1/$name1/filtered1/$name1.$b1.filtered1.vcf"
	fi
	if [ $start -le 12 ] && [ $end -ge 12 ]; then
		# clean stuff
		rm "$outdir1/$name1/$name1.2.bam"
		rm "$outdir1/$name1/$name1.3.bam"
		rm "$outdir1/$name1/$name1.4.bam"
		rm -r $outdir1"/"$name1"/TEMP"
	fi
fi

