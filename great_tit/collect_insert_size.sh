#!/bin/bash

set -u
set -e

if [ "$1" == -h ]; then
    echo ""
    echo ""
    echo ""
    echo ""
    echo ""
    echo "The first command line argument is a text file listing the full 
paths to the BAM files to be analysed on each row:"
    echo ""
    echo ""
    echo "example usage:"
    echo ""
    echo "sh collect_insert_size.sh bamfile_list.txt"
    echo ""
    exit 1
fi


#REF_FILE=/data/bo1pgc/gtit_1.04/genome/Parus_major_1.04.rename.fa # path to reference genome file that has been indexed 


#BGI files                
INSERT_METRICS=/data/bo1pgc/parus_reseq/bgi_reseq/insert_metrics #path where qc metrics will be written     
QSUB_PATH=/fastdata/bo1pgc/parus_reseq/bgi_reseq/qsub_out/insert_sizes_bgi

if [ ! -d "$INSERT_METRICS" ]; then
    mkdir $INSERT_METRICS
fi


if [ ! -d "$QSUB_PATH" ]; then
    mkdir $QSUB_PATH
fi


bam_files=$1

cut -f1 $bam_files | while read i;
do
    bam_id="$(basename $i .dedup.real.bam)"

    collect_inserts="java -Xmx2g -jar ~/bin/picard.jar CollectInsertSizeMetrics I=$i\
 O=$INSERT_METRICS/$bam_id.insert_metrics.txt HISTOGRAM_FILE=$INSERT_METRICS/$bam_id.insert_size.hist.pdf METRIC_ACCUMULATION_LEVEL=READ_GROUP"
    

    qsub -b y -N insert_metric.$bam_id -l mem=8G -l rmem=8G -l arch=intel* \
	-l h_rt=24:00:00 -e $QSUB_PATH/insert_metric.$bam_id.e \
	-o $QSUB_PATH/insert_metric.$bam_id.o -cwd $collect_inserts
done