#!/bin/bash

set -u
set -e


OUTDIR=/data/bo1pgc/parus_reseq/bgi_reseq/fastqc_report
QSUB_PATH=/fastdata/bo1pgc/parus_reseq/bgi_reseq/qsub_out/fastqc

if [ ! -d "$OUTDIR"  ]; then
    mkdir $OUTDIR
fi

if [ ! -d "$QSUB_PATH"  ]; then
    mkdir $QSUB_PATH
fi

queue='-q evolgen.q -P evolgen'


cut -f1,2 $1 | while read fq1 fq2; 
do
    bam_id="$(basename $fq1 .1.fq.gz )"
    
    qsub -b y $queue -N fastqc.$bam_id -l mem=8G -l rmem=8G -l arch=intel* \
        -l h_rt=8:00:00 -e $QSUB_PATH/fastqc.$bam_id.e \
        -o $QSUB_PATH/fastqc.$bam_id.o -cwd "~/bin/FastQC/fastqc $fq1 -o $OUTDIR"

    qsub -b y $queue -N fastqc.$bam_id -l mem=8G -l rmem=8G -l arch=intel* \
        -l h_rt=8:00:00 -e $QSUB_PATH/fastqc.$bam_id.e \
        -o $QSUB_PATH/fastqc.$bam_id.o -cwd "~/bin/FastQC/fastqc $fq2 -o $OUTDIR"

done


