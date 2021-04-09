
script_params='design/mapping_STAR.sh'
outdir="/Home/pH_experiment_adults/output/mapping_diamond"
db="/Home/prot5/merge/prot5" #selected proteomes

############################################

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

. "$script_params"
assert_

for i in ${!title[@]}; do
	echo "${title[i]}"
	echo "${for1[i]}"
	if [ -f "${for1[i]}" ]; then
		outIdx="$outdir/${title[i]}"
		q="${for1[i]}"
		echo "$outIdx"
		if [ $1 -eq 1 ]; then
		    # --more-sensitive
			if [ ! -f $outIdx.m8 ]; then
				echo OK
				t=20
				sbatch -N 1 -n"$t" --ntasks-per-node="$t" --mem=100000 -o "$outIdx.out" -e "$outIdx.err" \
				-p hive7d,hiveunlim,queen,ckpthp7d,ckptdell7d,ckptdell31d,ckpthp31d,preempt7d,ckptqueen,preempt31d,hive1d,ckpthp1d,preempt1d,ckptdell1d \
				--wrap "echo diamond
						function assert_ { rc=\$?; if [[ \$rc != 0 ]]; then echo 'exit !!!!'; exit \$rc; fi }
																
						diamond blastx -a \"$outIdx\" --tmpdir \"$outdir\" --threads \"$t\" -d \"$db\" -q \"$q\"  --index-chunks 1 --top 10 --evalue 0.01
						assert_
						diamond view -a $outIdx.daa -o $outIdx.m8; assert_
						rm $outIdx.daa; assert_
						
						"
			fi
		elif [ $1 -eq 2 ]; then
            if [ -f "$outIdx.m8" ]; then
                echo 'OK'
                sbatch -N1 -n1 --time=168:00:00 -p hive7d,ckpthp7d,ckpthp31d,preempt7d,preempt31d,ckptdell7d,ckptdell31d,hiveunlim,queen,ckptqueen \
                    --output="$outIdx.m8.krona.out" --error="$outIdx.m8.krona.err" \
                    --wrap "function assert_ { rc=\$?; if [[ \$rc != 0 ]]; then echo 'exit !!!!'; exit \$rc; fi }
                            cd \"$outdir\"; assert_
                            . \"$CONDA\"; assert_
                            conda activate krona_setting2; assert_ #2.7.1
                            blast1=\"$outIdx.m8\"
                            ktClassifyBLAST \"\$blast1\" -o \"\$blast1.krona.1.txt\"; assert_
                            ktImportBLAST   \"\$blast1.krona.1.txt\" -o \"\$blast1.krona.1.html\"; assert_
                            conda deactivate; assert_
                            "
            fi
		fi
	fi
	echo '----------------'
done



