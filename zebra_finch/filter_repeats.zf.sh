#!/bin/bash

#$ -l h_rt=08:00:00
#$ -l mem=8G
#$ -l rmem=8G
#$ -l arch=intel*
#$ -e filter_repeats.zf.e
#$ -o filter_repeats.zf.o
#$ -N filter_repeats
#$ -q evolgen.q -P evolgen


REPEATS=/data/bo1pgc/zebra_finch_genome/repeat_masking/zf_masking_coordinates.bed

java -Xmx3g -jar ~/bin/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R /data/bo1pgc/zebra_finch_genome/taeGut1.reordered.fa \
    -o $2 --variant $1 \
    --mask ${REPEATS} --maskName "REPEAT" 



