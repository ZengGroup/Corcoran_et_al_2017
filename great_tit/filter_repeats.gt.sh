#!/bin/bash

#$ -l h_rt=08:00:00
#$ -l mem=8G
#$ -l rmem=8G
#$ -l arch=intel*
#$ -e filter_repeats.gt.e
#$ -o filter_repeats.gt.o
#$ -N filter_repeats
#$ -q evolgen.q -P evolgen


REPEATS=/data/bo1pgc/gtit_1.04/repeat_masking/gt_masking_coordinates.bed

java -Xmx3g -jar ~/bin/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R ../ref_files/Parus_major_1.04.rename.fa \
    -o $2 --variant $1 \
    --mask ${REPEATS} --maskName "REPEAT" 



