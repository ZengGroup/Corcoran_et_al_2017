#/bin/bash

#This can be run in an interactive session

REF_PATH=/data/bo1pgc/gtit_1.04/genome

faReformat -s < $REF_PATH/Parus_major_1.04.fa > $REF_PATH/Parus_major_1.04.rename.fa # truncates records at first white space

# create the files needed downstream by GATK tools

java -Xmx2g -jar ~/bin/picard.jar CreateSequenceDictionary R=$REF_PATH/Parus_major_1.04.rename.fa O=$REF_PATH/Parus_major_1.04.rename.dict

bwa index $REF_PATH/Parus_major_1.04.rename.fa

samtools faidx $REF_PATH/Parus_major_1.04.rename.fa