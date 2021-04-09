#DE_file='/Home/STAR-spis-DE1/mesophotic/enrichment/mesophotic.comparison.txt'
#DE_file='/Home/STAR-spis-DE1/shallow/enrichment/shallow.comparison.txt'

#outDir='/Home/STAR-spis-DE1/shallow/enrichement/KEGG'
#geneToGo_file='/Home/kegg/spis.ptw.2.txt'
#title='KEGG'

#outDir='/Home/STAR-spis-DE1/shallow/enrichment/trinotate'
## geneToGo_file='/Home/trinotate/trinotate_annotation_report2.txt.geneLevel.GOs_.molten.2.withParents_'
#title='trinotate'

outDir='/Home/STAR-spis-DE1/shallow/enrichment/uniprot'
geneToGo_file='/Home/uniprot/uniprot-stylophora.tab.GOs.molten_gene-to-GOterms1.txt'
title='uniprot'

column_id='ID'
start=1
end=20

######################################

function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }

if [ ! -d "$outDir" ] || [ ! -f "$DE_file" ] ||  [ ! -f "$geneToGo_file" ]; then echo "not all files or directories"; exit; fi
if [ ! -d "$outDir/log" ]; then mkdir "$outDir/log"; assert_; fi

if [ $1 -eq 0 ]; then
    head -n1 "$DE_file"
elif [ $1 -eq 1 ]; then
    for i in $(seq $start $end); do
        echo "i = $i"
        sbatch -o "$outDir/log/$title.$i.out" -e "$outDir/log/$title.$i.err" \
            -p hive1d,hive7d,ckpthp1d,ckpthp7d,ckpthp31d,preempt1d,preempt7d,preempt31d,ckptdell1d,ckptdell7d,ckptdell31d,hiveunlim,queen,ckptqueen \
            --wrap ". \"$CONDA\"
			conda activate R_setting6
			Rscript goseq_non-natively-supported-108.R \
            --start=$i \
            --end=$i \
             --outDir=\"$outDir\" \
            --DE_file=\"$DE_file\" \
            --geneToGo_file=\"$geneToGo_file\" \
            --title=\"$title\" \
            --column_id=\"$column_id\" \
            --splitID=\"F\"
			conda deactivate"
		sleep 2
        #break
    done
elif [ $1 -eq 2 ]; then
    
    cd "$outDir"
    summary_file_1_=("$title.goseq_nobias.summary.txt" "$title.goseq_wall.summary.txt")
    search_=("goseq_nobias.txt" "goseq_wall.txt")
    
    for i in "${!summary_file_1_[@]}"; do
        summary_file_1="${summary_file_1_[$i]}"
        search1="${search_[$i]}"
        echo "writing into \"$summary_file_1\" ..."
        if [ -f "$summary_file_1" ]; then rm "$summary_file_1"; fi
        visited=0
        for z in */"$title"*."$search1"; do
            n=$(wc -l "$z" | sed 's/ .*//')
            if [ $n -gt 1 ]; then
                name1=$(echo "$z" | sed 's/.*\///' | sed 's/.goseq_.*//' | perl -ne '$_ =~ s/\.select\.\d+$|^\w+\.//g; printf("%s",$_)')
                echo "$name1"
                if [ $visited -eq 0 ]; then
                    cat "$z" | \
                        perl -slane 'chomp; if($. > 1){printf("%s\t%s\n",$name1,$_)}else{printf("%s\t%s\n","comparison",$_)}' -- -name1="$name1" | \
                        perl -ne 'chomp; @x=split /\t/,$_; for($i=0;$i<$#x+1;$i++){$x[$i]="\"".$x[$i]."\""} printf("%s\n",join("\t",@x))' \
                        >> "$summary_file_1"
                else
                    cat "$z" | \
                        perl -slane 'chomp; if($. > 1){printf("%s\t%s\n",$name1,$_)}' -- -name1="$name1" | \
                        perl -ne 'chomp; @x=split /\t/,$_; for($i=0;$i<$#x+1;$i++){$x[$i]="\"".$x[$i]."\""} printf("%s\n",join("\t",@x))' \
                        >> "$summary_file_1"
                fi
                
                visited=1
            fi
        done
    done
	
    zip -j "$title.goseq_fit.zip" */"$title"*.pdf

fi

