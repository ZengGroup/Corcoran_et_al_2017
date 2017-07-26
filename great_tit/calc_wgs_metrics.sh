#!/bin/bash


REF_FILE=/data/bo1pgc/gtit_1.04/genome/Parus_major_1.04.rename.fa # path to reference genome file that has been indexed
BAM_FILES=/fastdata/bo1pgc/parus_reseq/bgi_reseq/recalibrated_bams # path where recalibrated bam files are located             
QC_METRICS=/data/bo1pgc/parus_reseq/bgi_reseq/mapping_qc_recal_bams #path where qc metrics will be written           
QSUB_PATH=/fastdata/bo1pgc/parus_reseq/bgi_reseq/qsub_out/wgs_metric_recal_bams # directory to write the qsub job log files                                                 

if [ ! -f "$REF_FILE" ];  then
    echo "The reference genome file path is incorrect"
    exit 1
fi

if [ ! -d "$BAM_FILES" ]; then
    echo "The path to the BAM files is incorrect"
    exit 1
fi

if [ ! -d "$QC_METRICS" ]; then
    mkdir $QC_METRICS
fi


if [ ! -d "$QSUB_PATH" ]; then
    mkdir $QSUB_PATH
fi

                                                                                                   
bam_files=$( ls $BAM_FILES/*.bam )

for i in $bam_files;
do
    bam_id="$(basename $i .dedup.real.recal.bam)"
    
    
    wgs_metrics="java -Xmx2g -jar ~/bin/picard.jar CollectWgsMetrics I=$i 
O=$QC_METRICS/$bam_id.wgsmetrics_file.txt R=$REF_FILE INCLUDE_BQ_HISTOGRAM=true "
    
    qsub -b y -N wgs_metrics.$bam_id -l mem=8G -l rmem=8G -l arch=intel* \
        -l h_rt=24:00:00 -e $QSUB_PATH/wgs_metrics.$bam_id.e \
        -o $QSUB_PATH/wgs_metrics.$bam_id.o -cwd $wgs_metrics

done
