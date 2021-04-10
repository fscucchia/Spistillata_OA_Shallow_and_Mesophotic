
path1=/Home/STAR-spis/*.starAligned.sortedByCoord.out.bam
fastqPhred=33
remove_from_title=".starAligned.sortedByCoord.out.bam"
gtf1='/Home/data/stylophora-genome/liftover3-NCBI-spis/merged-NCBI-GCF_002571385.1_and_spis.sorted.FIX1_.gtf'

g1='/Home/data/stylophora-genome/GCF_002571385.1_Stylophora_pistillata_v1_genomic.fna'
outdir1='/Home/STAR-spis-SNPs/gatk1'

GVCF_outDir1='/Home/STAR-spis-SNPs/gatk1/joint_analysis'
gvcfs_title='pH2'

whichGenotyper=1 # HaplotypeCaller & HaplotypeCaller GVCF = 1, test = 0
bamVersion=5 # running haplotype caller with bam file version: 5 (GATK original setting), or 6 (after removing all multireads with bamutils)
annotationDB="Stylophora_pistillata_v1_GCF_002571385" # snpEff database
vcfRefDB="" # "" or vcf
ploidy=2

#java7exe=/data/apps/java/jdk1.7.0_79_64bit/bin/java
java7exe=/Home/programs/jdk1.8.0_65/bin/java
gatkExe=/Home/programs/gatk/GATK3.5 # the directory where GenomeAnalysisTK.jar is located
snpEff=/Home/programs/snpEff_ # the directory where snpEff.jar is in
picard=/Home/programs/picard-tools-1.140/picard.jar

# VariantFiltration params:
#cond1="DP < 45 || QUAL < 150.0"
cond1="DP < 6.0"
cond1mame="DP6"
#min_af=0.2
#max_af=0.8

################################################################

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

if [ $# -eq 3 ]; then
	# 3 input arguments
	part=$1      # part=1 for the GATK pipline, part=2 for GATK output processing
	stepStart=$2 # start step
	stepEnd=$3   # end step
else
	echo "incorrect numbers of args !!!"
	exit 1
fi

i=0
if [ $part -eq 1 ] || [ $part -eq 2 ]; then
	for a1 in $path1; do 
		echo $a1
		#name1=`echo $a1 | grep -P -o $name_regexpr`
        name1=$(echo $a1 | sed 's/.*\///' | sed 's/'$remove_from_title'//')
		echo $name1
		let i=$i+1
		echo $i
	
		if [ $part -eq 1 ]; then
			if [ $stepStart -eq 0 ]; then
				# creating a fasta index and dict for the reference genome, this is done once
				sbatch -N1 -n1 --ntasks-per-node=1 \
                    -o "$outdir1/$name1.out" -e "$outdir1/$name1.err"  \
                    -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d \
                    --wrap "bash gatk21_pipeline.sh \"preserved\" 0 0 \
							\"$a1\" \"$g1\" \"$outdir1\" \"$whichGenotyper\" \"$name1\" \
							\"$fastqPhred\" \"$bamVersion\" \"$ploidy\" \"$vcfRefDB\" \
							\"$picard\" \"$gatkExe\" \"$java7exe\" \"$snpEff\" 1 1"
                    
				break # don't remove the break command here
			elif [ $stepStart -ge 1 ] && [ $stepEnd -le 12 ]; then
				#if [ ! -d $outdir1/$name1 ] || [ $stepStart -eq 11 ] || [ $stepStart -eq 12 ]; then
					echo "$i) $stepStart"
					echo "$i) $stepEnd"
					sbatch -N1 -n1 --ntasks-per-node=1 --mem=64000 --time=10080 \
                        -o "$outdir1/$name1.out" -e "$outdir1/$name1.err"  \
                        -p hive7d,hiveunlim,queen,preempt7d,preempt31d,ckpthp7d,ckpthp31d,ckptdell7d,ckptdell31d \
                            --wrap "bash gatk21_pipeline.sh \"preserved\" $stepStart $stepEnd \
                                    \"$a1\" \"$g1\" \"$outdir1\" \"$whichGenotyper\" \"$name1\" \
                                    \"$fastqPhred\" \"$bamVersion\" \"$ploidy\" \"$vcfRefDB\" \
                                    \"$picard\" \"$gatkExe\" \"$java7exe\" \"$snpEff\" 1 1"
				#fi
			fi
		elif [ $part -eq 2 ]; then
			echo $outdir1/$name1/$name1.5.g.vcf
			echo $outdir1/$name1/$name1.5.vcf
			echo $outdir1/$name1/filtered1/$name1.5.filtered1.vcf
			if [ -e $outdir1/$name1/$name1.5.g.vcf ] && [ -e $outdir1/$name1/$name1.5.vcf ] && [ -e $outdir1/$name1/filtered1/$name1.5.filtered1.vcf ]; then
				echo "--> $part $stepStart $stepEnd"
				vcf1=$outdir1/$name1/filtered1/$name1.5.filtered1.vcf
				echo $vcf1
                
                sbatch -N1 -n1 --ntasks-per-node=1 --mem=64000 \
                        -o "$outdir1/$name1.out" -e "$outdir1/$name1.err"  \
                        -p hive7d,hiveunlim,queen,preempt7d,preempt31d,ckpthp7d,ckpthp31d,ckptdell7d,ckptdell31d \
				        --wrap "bash gatk21_vcf-processing.sh \"preserved\" \"$stepStart\" \"$stepEnd\" \
									\"$vcf1\" \"$g1\" \"$annotationDB\" \"$cond1\" \"$cond1mame\" \
									\"$vcfRefDB\" \"$picard\" \"$gatkExe\" \"$java7exe\" \"$snpEff\" \"$gtf1\""
               # \"$min_af\" \"$max_af\" 
			fi
		fi
	
		echo '----------------------------'

	done
fi
#elif [ $part -eq 3 ]; then
#	if [ $stepStart -ge 1 ] && [ $stepEnd -le 1 ]; then
#		#
#		# writes a script
#		GVCFs="$outdir1/*/filtered1/*.snps.pass.vcf"
#        echo "$GVCFs"
#		python GenotypeGVCFs5.py "$java7exe" "$gatkExe" "$g1" "$GVCFs" "$GVCF_outDir1" "$ploidy" 20 "$gvcfs_title"
#		assert_
#        echo ""
#		echo "check and run:"
#		echo "$GVCF_outDir1/$gvcfs_title.joint_analysis.sh"
#		
#	elif [ $stepStart -ge 2 ] && [ $stepEnd -le 2 ] && [ -e "$GVCF_outDir1/$gvcfs_title.joint_analysis.sh" ]; then
#		
#		echo "$GVCF_outDir1/$gvcfs_title.joint_analysis.sh"
#		
#		sbatch SBATCH_p20_ckpt7d_mem128_1.sh "$GVCF_outDir1/$gvcfs_title.joint_analysis.sh" 1 &
#		#sbatch SBATCH_p32_ckpt7d_mem764_1.sh "$GVCF_outDir1/$gvcfs_title.joint_analysis.sh" 1 &
#	
#	elif [ $stepStart -ge 3 ] && [ $stepEnd -le 3 ] && [ -e "$GVCF_outDir1/joint_data_$gvcfs_title.vcf" ]; then
#		#z='.ncbi-annotations.test1'
#		z=''
#		echo "$GVCF_outDir1/joint_data_$gvcfs_title.vcf"
#		if [ ! -e "$GVCF_outDir1/filtered1$z/joint_data_$gvcfs_title.vcf" ]; then
#			mkdir "$GVCF_outDir1/filtered1$z"
#			cp "$GVCF_outDir1/joint_data_$gvcfs_title.vcf" "$GVCF_outDir1/filtered1$z/joint_data_$gvcfs_title.vcf"
#			assert_
#		fi
#		echo "$GVCF_outDir1/filtered1$z/joint_data_$gvcfs_title.vcf"
#		vcf1="$GVCF_outDir1/filtered1$z/joint_data_$gvcfs_title.vcf"
#	fi
#fi

