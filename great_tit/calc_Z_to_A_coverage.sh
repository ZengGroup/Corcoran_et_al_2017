#!/bin/bash

set -e
set -u

BAM_FILES=$1
GENOME_FILE=$2
OUTFILE=$3


QSUB_PATH=/fastdata/bo1pgc/parus_reseq/bgi_reseq/qsub_out/Z_to_aut_coverage

if [ ! -d "$QSUB_PATH" ]; then
    mkdir $QSUB_PATH
fi

if [ -f "$OUTFILE" ]; then
    rm $OUTFILE
fi

cat $BAM_FILES | while read i j
do
    
    sample_id="$(basename $i .dedup.real.bam)"
    
    samtools_depth="samtools view -h -u -F 1024 $i | samtools depth                                                                          
 -q20 -Q20 - | python calc_Z_to_A_coverage.py $i $2 $3  >> $OUTFILE"

    qsub -b y -N calc_depth.$sample_id -l mem=4G -l rmem=4G  -l arch=intel* \
       -l h_rt=8:00:00 -e $QSUB_PATH/calc_depth.$sample_id..e \
       -o $QSUB_PATH/calc_depth.$sample_id.o -cwd $samtools_depth

done