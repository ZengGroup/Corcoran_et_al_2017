#!/bin/bash

set -u
set -e

if [ "$1" == -h ]; then
    echo ""
    echo "This script call samtools flagstat command to get a summary of 
the percentage reads mapped and unmapped"
    echo ""
    echo ""
    echo ""
    echo "The first command line argument is a text file listing the full 
paths to the BAM files to be analysed on each row:"
    echo ""
    echo ""
    echo "Example usage:"
    echo ""
    echo "sh get_mapping_stats.sh bamfile_list.txt"
    echo ""
    exit 1
fi



MAPPING_METRICS=/data/bo1pgc/parus_reseq/bgi_reseq/samtools_flagstat #path where mapping metrics will be written     
QSUB_PATH=/fastdata/bo1pgc/parus_reseq/bgi_reseq/qsub_out/samtools_flagstat

if [ ! -d "$MAPPING_METRICS" ]; then
    mkdir $MAPPING_METRICS
fi


if [ ! -d "$QSUB_PATH" ]; then
    mkdir $QSUB_PATH
fi


bam_files=$1

cut -f1 $bam_files | while read i;
do
    bam_id="$(basename $i .dedup.real.recal.bam)"
    
    samtools_flagstat="samtools flagstat $i > $MAPPING_METRICS/$bam_id.flagstat.txt"

    qsub -b y -N mapping_metric.$bam_id -l mem=6G -l rmem=6G -l arch=intel* \
	-l h_rt=8:00:00 -e $QSUB_PATH/mapping_metric.$bam_id.e \
	-o $QSUB_PATH/mapping_metric.$bam_id.o -cwd $samtools_flagstat

done