#!/bin/bash
#$ -o /u/scratch/a/asperlea/joblogs/
#$ -e /u/scratch/a/asperlea/joblogs/
#$ -m a
#$ -pe shared 12
#$ -l h_data=4G,highp,h_rt=168:00:00

. /u/local/Modules/default/init/modules.sh
module load java
java -jar /u/project/ernst/asperlea/ConservationProject/ChromHMM/ChromHMM.jar LearnModel -b 1 -nobed -n 150 -d -1 -lowmem -p 12 /u/project/ernst/asperlea/ConsHMM_annotations/binarizedFiles/hg38_multiz100way /u/project/ernst/asperlea/ConsHMM_annotations/models/hg38_multiz100way/100_states/ 100 hg19
