
bed1=/data/home/STAR-spis-SNPs/gatk1-joint/pH2-joint-all.2.bed
outdir=/data/home//STAR-spis-SNPs/gatk1-joint/mosDepth
bams_path=/data/home//STAR-spis-SNPs/gatk1/*/*.5.bam

if [ $1 -le 2 ]; then
    for bam1 in $bams_path; do
        echo $bam1
        x=$(echo "$bam1" | sed 's/.*\///' | sed 's/.5.bam//')
        echo "$x"
        prefix1="$outdir"/"$x"
        echo $prefix1
        
        if [ -f "$bed1" ] & [ -d "$outdir" ]; then
            echo 'OK'
            if [ $1 -eq 2 ]; then
                sbatch --mem=128000 -p hiveunlim,hive1d,hive7d,queen -N1 -n4 \
                    -o "$prefix1.out" -e "$prefix1.err" \
                    --wrap "
                            if [ ! -z ${CONDA+x} ]; then
                                function assert_ { rc=\$?; if [[ \$rc != 0 ]]; then echo 'exit !!!!'; exit \$rc; fi }
                                . \"$CONDA\"; assert_
                                #conda activate mosdepth_1; assert_
                                conda activate mosdepth_setting1; assert_
                                echo 'cnda is activated'
                                mosdepth --fast-mode --no-per-base --threads 4 --by \"$bed1\" \"$prefix1\" \"$bam1\"; assert_
                                conda deactivate; assert_
                            else
                                echo \"CONDA global variable doesn't exist\"; exit
                            fi
                        "
            fi
        fi
        echo '--------'
        #break
    done
elif [ $1 -eq 3 ]; then
    gunzip "$outdir"/*.regions.bed.gz
elif [ $1 -eq 4 ]; then
    for x in "$outdir"/*.regions.bed; do
        echo "$x"
        y=$(echo "$x" | sed 's/.*\///' | sed 's/.regions.bed//')
        echo "$y"
        cat "$x" | awk "BEGIN{FS=\" \"; print \"ID.$y coverage.$y\"} {print \$4,\$5}" > "$x.1"
        echo '---------'
    done
elif [ $1 -eq 5 ]; then
    echo 'xxx'
    sbatch --mem=128000 -p hiveunlim,hive1d,hive7d,queen -N1 -n1 -o "$outdir/merge.out" \
        --wrap "ls \"$outdir\"/*.regions.bed.1 | xargs paste -d ' ' > \"$outdir\"/merge.txt"
fi

