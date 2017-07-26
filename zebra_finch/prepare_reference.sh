#!/usr/bin/env bash

#$ -l h_rt=01:00:00
#$ -l mem=4G
#$ -l rmem=4G
#$ -l arch=intel*
#$ -e /fastdata/bo1pgc/parus_reseq/singhal_zebra_finch_data/zf_snp_calling/qsub_out/prepare_ref.e
#$ -o /fastdata/bo1pgc/parus_reseq/singhal_zebra_finch_data/zf_snp_calling/qsub_out/prepare_ref.o
#$ -N prepare_ref


#This can be run in an interactive session

REF_PATH=/data/bo1pgc/zebra_finch_genome/


# create the files needed downstream by GATK tools

java -Xmx2g -jar ~/bin/picard.jar CreateSequenceDictionary R=$REF_PATH/taeGut1.reordered.fa  O=$REF_PATH/taeGut1.reordered.dict

bwa index $REF_PATH/taeGut1.reordered.fa 

samtools faidx $REF_PATH/taeGut1.reordered.fa 