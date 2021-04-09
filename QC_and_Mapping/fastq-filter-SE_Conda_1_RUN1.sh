
design_file='design/design1.concatenated.sh'
dir2='$HOME1/pH_experiment_adults/output/filtered_raw'

##################################################################

# check jobs status: sacct --state r ; sacct --state cd
# run on 5 nodes, totally 50 threads, 10 threads per node: srun -N5 -n50 --ntasks-per-node=10

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

mkdir -p $dir2"/fastqc"; assert_
mkdir -p $dir2"/fastq_screen"; assert_

. "$design_file"; assert_

if [ $1 -le 4 ]; then
	for i in ${!title[@]}; do
		f1="${for1[i]}"
		title="${title[i]}"
		f1_=$dir2"/"$title".filtered"

		#if [ ! -f $f1_ ]; then
		
		echo '----------------'
		echo $i
		echo ${title[i]}
		echo ${for1[i]}
		echo $f1_
		
		if [ $1 -eq 1 ]; then
			echo '1'
			sbatch -N1 -n1 --ntasks-per-node=1 --workdir=$dir2/fastqc \
				-p hive1d,hive7d,ckpthp1d,ckpthp7d,ckpthp31d,preempt1d,preempt7d,preempt31d,ckptdell1d,ckptdell7d,ckptdell31d,hiveunlim,queen,ckptqueen \
				--wrap ". $CONDA; conda activate rnaseq_1; fastqc -o $dir2/fastqc $f1; conda deactivate"
		elif [ $1 -eq 2 ]; then
			echo '2'
			qual_threshold1=25
			qual_window1=10
			minlen=20
			threads=20
			sbatch -N1 -n20 --ntasks-per-node=20 --mem=128000 -o $dir2/$title.out -e $dir2/$title.err \
				-p hive1d,hive7d,ckpthp1d,ckpthp7d,ckpthp31d,preempt1d,preempt7d,preempt31d,ckptdell1d,ckptdell7d,ckptdell31d,hiveunlim,queen,ckptqueen \
				fastq-filter_Conda_job_1.sh \
					'SE' 'adapters4d.fa' $qual_threshold1 $qual_window1 $minlen $threads \
					$f1 $f1_
		elif [ $1 -eq 3 ]; then
			echo '3'
			sbatch -N1 -n1 --ntasks-per-node=1 --workdir=$dir2/fastqc -o $dir2/$title.out -e $dir2/$title.err \
				-p hive1d,hive7d,ckpthp1d,ckpthp7d,ckpthp31d,preempt1d,preempt7d,preempt31d,ckptdell1d,ckptdell7d,ckptdell31d,hiveunlim,queen,ckptqueen \
				--wrap ". $CONDA; conda activate rnaseq_1; fastqc -o $dir2/fastqc $f1_; conda deactivate"
		elif [ $1 -eq 4 ]; then
			echo '4'
			if [ -d $dir2"/fastq_screen" ]; then
				sbatch -N1 -n20 --ntasks-per-node=20 --mem=128000 \
					-e $dir2"/fastq_screen"/$title.err -o $dir2"/fastq_screen"/$title.out \
					-p hive1d,hive7d,ckpthp1d,ckptdell1d,ckpthp7d,ckpthp31d,preempt7d,preempt31d,ckptdell7d,ckptdell31d,queen \
						fastq_screen --nohits --outdir $dir2"/fastq_screen" --conf fastq_screen.conf \
										--aligner bowtie2 --threads 20 $f1_
			fi
		fi
		
		#break
		#fi
	done
elif [ $1 -eq 5 ]; then
	echo '5'
	sbatch -N1 -n1 --ntasks-per-node=1 --workdir=$dir2/fastqc -o "$dir2/fastqc/multiqc.out" -e  "$dir2/fastqc/multiqc.err" \
		-p hive1d,hive7d,ckpthp1d,ckpthp7d,ckpthp31d,preempt1d,preempt7d,preempt31d,ckptdell1d,ckptdell7d,ckptdell31d,hiveunlim,queen,ckptqueen \
		--wrap ". $CONDA; conda activate rnaseq_1; cd \"$dir2/fastqc\"; multiqc .; conda deactivate"
fi
