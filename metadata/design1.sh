# comments:
# script template created by mapper_write-SE-or-PE-list_12.py from table design/design1.csv

# files:
function assert_ { rc=$?; if [[ $rc != 0 ]]; then echo 'exit !!!!'; exit $rc; fi }
title[0]="A1"
for1[0]="/Home/fastq/A1/A1_S1_L002_R1_001.fastq.gz,/Home/fastq/A1/A1_S1_L001_R1_001.fastq.gz"
title[1]="A4"
for1[1]="/Home/fastq/A4/A4_S2_L001_R1_001.fastq.gz,/Home/fastq/A4/A4_S2_L002_R1_001.fastq.gz"
title[2]="A7"
for1[2]="/Home/fastq/A7/A7_S3_L002_R1_001.fastq.gz,/Home/fastq/A7/A7_S3_L001_R1_001.fastq.gz"
title[3]="B1"
for1[3]="/Home/fastq/B1/B1_S4_L001_R1_001.fastq.gz,/Home/fastq/B1/B1_S4_L002_R1_001.fastq.gz"
title[4]="B4"
for1[4]="/Home/fastq/B4/B4_S5_L001_R1_001.fastq.gz,/Home/fastq/B4/B4_S5_L002_R1_001.fastq.gz"
title[5]="B7"
for1[5]="/Home/fastq/B7/B7_S6_L002_R1_001.fastq.gz,/Home/fastq/B7/B7_S6_L001_R1_001.fastq.gz"
title[6]="C3"
for1[6]="/Home/fastq/C3/C3_S7_L001_R1_001.fastq.gz,/Home/fastq/C3/C3_S7_L002_R1_001.fastq.gz"
title[7]="C4"
for1[7]="/Home/fastq/C4/C4_S8_L001_R1_001.fastq.gz,/Home/fastq/C4/C4_S8_L002_R1_001.fastq.gz"
title[8]="C7"
for1[8]="/Home/fastq/C7/C7_S9_L001_R1_001.fastq.gz,/Home/fastq/C7/C7_S9_L002_R1_001.fastq.gz"
title[9]="V2_1"
for1[9]="/Home/fastq/V2_1/V2_1_S34_L002_R1_001.fastq.gz,/Home/fastq/V2_1/V2_1_S34_L001_R1_001.fastq.gz"
title[10]="VI_2_1"
for1[10]="/Home/fastq/VI_2_1/VI_2_1_S35_L002_R1_001.fastq.gz,/Home/fastq/VI_2_1/VI_2_1_S35_L001_R1_001.fastq.gz"
title[11]="VI_2_4"
for1[11]="/Home/fastq/VI_2_4/VI_2_4_S36_L002_R1_001.fastq.gz,/Home/fastq/VI_2_4/VI_2_4_S36_L001_R1_001.fastq.gz"
title[12]="VI_2_7"
for1[12]="/Home/fastq/VI_2_7/VI_2_7_S37_L001_R1_001.fastq.gz,/Home/fastq/VI_2_7/VI_2_7_S37_L002_R1_001.fastq.gz"
title[13]="V_2_4"
for1[13]="/Home/fastq/V_2_4/V_2_4_S38_L001_R1_001.fastq.gz,/Home/fastq/V_2_4/V_2_4_S38_L002_R1_001.fastq.gz"
title[14]="V_2_7"
for1[14]="/Home/fastq/V_2_7/V_2_7_S39_L002_R1_001.fastq.gz,/Home/fastq/V_2_7/V_2_7_S39_L001_R1_001.fastq.gz"
title[15]="X_2_1"
for1[15]="/Home/fastq/X_2_1/X_2_1_S40_L002_R1_001.fastq.gz,/Home/fastq/X_2_1/X_2_1_S40_L001_R1_001.fastq.gz"
title[16]="X_2_4"
for1[16]="/Home/fastq/X_2_4/X_2_4_S41_L002_R1_001.fastq.gz,/Home/fastq/X_2_4/X_2_4_S41_L001_R1_001.fastq.gz"
title[17]="X_2_8"
for1[17]="/Home/fastq/X_2_8/X_2_8_S42_L001_R1_001.fastq.gz,/Home/fastq/X_2_8/X_2_8_S42_L002_R1_001.fastq.gz"

#run:
#for i in ${!title[@]}; do
#	echo $i
#	echo ${title[i]}
#	echo ${for1[i]}
#	echo '----------------'
#done
