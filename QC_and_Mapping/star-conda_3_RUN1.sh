
script_params='design/mapping_STAR.sh'
genome1='$HOME1/$WD1/data/stylophoraP.fasta'
star_genome_dir='$HOME1/$WD1/data/STAR_spis-liftover3'
gtf='$HOME1/$WD1/data/dbs/stylophoraP-liftover3.gtf'
gtfGeneReg="(?<!\")smic_gene\d+(?=\")|(?<!\")gene\d+(?=\")|(?<!\")SpisGene\d+(?=\")"
outdir='$HOME1/pH_experiment_adults/output/mapping'
pairing="SE"
# min(14, log2(GenomeLength)/2 - 1) # 400megabase genome size gives 13, default  14
genomeSAindexNbases=14 # 808226712 bp, 13.7950924
limitSjdbInsertNsj=1200000 # default 1000000

#gffread "$gff1" -E -T -F -o "$gff1.gtf" | tee  "$gff1.gtf.errors"
#gffread "$gff1" -g "$genome1" -E -F -w "$gff1.rna.fasta"
#gffread "$gff1" -g "$genome1" -E -F -x "$gff1.CDs.fasta"

###########################################

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. $script_params
assert_

if [ ! -z ${CONDA+x} ]; then
	if [ $1 -eq 1 ]; then
		mkdir -p "$star_genome_dir"; assert_

		sbatch --mem=128000 -N1 -n20 --ntasks-per-node=20 \
		    -o "$star_genome_dir/star.out" -e "$star_genome_dir/star.err" \
		    -p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d \
		    --wrap "
					if [ ! -z ${CONDA+x} ]; then
						function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }
						. \"$CONDA\"; assert_
						conda activate star2.5; assert_
						STAR --runMode genomeGenerate --runThreadN 20 --genomeFastaFiles \"$genome1\" --genomeDir \"$star_genome_dir\" \
							--sjdbGTFfile \"$gtf\" --outFileNamePrefix \"$star_genome_dir/star\" --genomeSAindexNbases $genomeSAindexNbases \
							--limitSjdbInsertNsj $limitSjdbInsertNsj
						assert_
						conda deactivate; assert_
					else
						echo \"CONDA global variable doesn't exist\"; assert_
					fi
					"
	elif [ $1 -ge 2 ] && [ $1 -le 3 ] && [ -d "$outdir" ] && [ -f "$gtf" ]; then
		for i in ${!title[@]}; do
            t="$outdir/${title[i]}"
			if [ ! -f "$t.starReadsPerGene.out.tab" ]; then
				echo "title = ""${title[i]}"
				echo "for = ""${for1[i]}"
				if [ ! -z "${rev1[i]}" ]; then echo "rev = ${rev1[i]}"; fi
				if [ $1 -eq 3 ]; then
					if [ $pairing == "SE" ]; then
						if [ ! -z "${for1[i]}" ] && [ ! -z "${title[i]}" ]; then
							if [ -f "${for1[i]}" ]; then
								sbatch --mem=128000 -N1 -n20 --ntasks-per-node=20 \
									-o "$t".star.out -e "$t".star.err \
									-p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d \
									--wrap "echo SE
											function assert_ { rc=\$?; if [[ \$rc != 0 ]]; then echo 'exit !!!!'; exit \$rc; fi }
											if [ ! -z ${CONDA+x} ]; then
												. \"$CONDA\"; assert_
												conda activate star2.5; assert_
												STAR --readFilesIn \"${for1[i]}\" --outFileNamePrefix \"$t.star\" --runThreadN 20 --genomeDir \"$star_genome_dir\" \
													--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
												assert_
												conda deactivate; assert_
											else
												echo \"CONDA global variable doesn't exist\"; assert_
											fi
											"
							fi
						fi
					elif [ $pairing == "PE" ]; then
						if [ ! -z "${for1[i]}" ] && [ ! -z "${rev1[i]}" ] && [ ! -z "${title[i]}" ]; then
							if [ -f "${for1[i]}" ] && [ -f "${rev1[i]}" ]; then
								sbatch --mem=128000 -N1 -n20 --ntasks-per-node=20 \
									-o "$t".star.out -e "$t".star.err \
									-p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d \
									--wrap "echo PE
											function assert_ { rc=\$?; if [[ \$rc != 0 ]]; then echo 'exit !!!!'; exit \$rc; fi }
											if [ ! -z ${CONDA+x} ]; then
												. \"$CONDA\"; assert_
												conda activate star2.5; assert_
												STAR --readFilesIn \"${for1[i]}\" \"${rev1[i]}\" --outFileNamePrefix \"$t.star\" --runThreadN 20 --genomeDir \"$star_genome_dir\" \
													--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
												assert_
												conda deactivate; assert_
											else
												echo \"CONDA global variable doesn't exist\"; assert_
											fi
											"
							fi
						fi
					fi
				fi
				echo "----------------"
			fi
		    #break
		done
	elif [ $1 -eq 4 ]; then
		python gene-info-from-gtf_1.py --gtf "$gtf" --gene_regexpr "$gtfGeneReg" --output_file "$gtf.gene-length.txt"
	elif [ $1 -eq 5 ]; then
		mkdir -p "$outdir/multiqc"; assert_

		sbatch --mem=128000 -N1 -n20 --ntasks-per-node=20 \
			-o "$outdir/multiqc/multiqc".out -e "$outdir/multiqc/multiqc".err \
			-p hive1d,hive7d,hiveunlim,queen,preempt1d,preempt7d,preempt31d,ckpthp1d,ckpthp7d,ckpthp31d,ckptdell1d,ckptdell7d,ckptdell31d \
				--wrap "
						function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }
						if [ ! -z ${CONDA+x} ]; then
							. \"$CONDA\"; assert_
							conda activate rnaseq_1; assert_
							cd $outdir/multiqc
							multiqc $outdir
							assert_
							conda deactivate; assert_
						else
							echo \"CONDA global variable doesn't exist\"; assert_
						fi
						"
	fi
else
	echo "CONDA global variable doesn't exist"; assert_
fi

