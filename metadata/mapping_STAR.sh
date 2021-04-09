
# files:
function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }
title[0]="A1"
for1[0]="/Home/pH_experiment_adults/output/filtered_raw/A1.filtered"
title[1]="A4"
for1[1]="/Home/pH_experiment_adults/output/filtered_raw/A4.filtered"
title[2]="A7"
for1[2]="/Home/pH_experiment_adults/output/filtered_raw/A7.filtered"
title[3]="B1"
for1[3]="/Home/pH_experiment_adults/output/filtered_raw/B1.filtered"
title[4]="B4"
for1[4]="/Home/pH_experiment_adults/output/filtered_raw/B4.filtered"
title[5]="B7"
for1[5]="/Home/pH_experiment_adults/output/filtered_raw/B7.filtered"
title[6]="C3"
for1[6]="/Home/pH_experiment_adults/output/filtered_raw/C3.filtered"
title[7]="C4"
for1[7]="/Home/pH_experiment_adults/output/filtered_raw/C4.filtered"
title[8]="C7"
for1[8]="/Home/pH_experiment_adults/output/filtered_raw/C7.filtered"
title[9]="V2_1"
for1[9]="/Home/pH_experiment_adults/output/filtered_raw/V2_1.filtered"
title[10]="VI_2_1"
for1[10]="/Home/pH_experiment_adults/output/filtered_raw/VI_2_1.filtered"
title[11]="VI_2_4"
for1[11]="/Home/pH_experiment_adults/output/filtered_raw/VI_2_4.filtered"
title[12]="VI_2_7"
for1[12]="/Home/pH_experiment_adults/output/filtered_raw/VI_2_7.filtered"
title[13]="V_2_4"
for1[13]="/Home/pH_experiment_adults/output/filtered_raw/V_2_4.filtered"
title[14]="V_2_7"
for1[14]="/Home/pH_experiment_adults/output/filtered_raw/V_2_7.filtered"
title[15]="X_2_1"
for1[15]="/Home/pH_experiment_adults/output/filtered_raw/X_2_1.filtered"
title[16]="X_2_4"
for1[16]="/Home/pH_experiment_adults/output/filtered_raw/X_2_4.filtered"
title[17]="X_2_8"
for1[17]="/Home/pH_experiment_adults/output/filtered_raw/X_2_8.filtered"

#run:
#for i in ${!title[@]}; do
#	echo $i
#	echo ${title[i]}
#	echo ${for1[i]}
#	echo '----------------'
#done
