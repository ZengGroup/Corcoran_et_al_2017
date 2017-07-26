#!/bin/bash

QSUB_PATH=/fastdata/bo1pgc/parus_reseq/singhal_zebra_finch_data/zf_snp_calling/qsub_out

bam_list=bamlist.list
bam_array=($(cat $bam_list)) 


for i in ${bam_array[@]}
do
    sample_id="$(basename $i .recal.bam)"
    
    index_command="java -Xmx4g -jar ~/bin/picard.jar BuildBamIndex INPUT=$i"

    qsub -b y -N index_bam.$sample_id -l mem=6G -l rmem=6G  -l arch=intel* \
	    -l h_rt=1:00:00 -e $QSUB_PATH/bam_index.$sample_id.e \
	    -o $QSUB_PATH/bam_index.$sample_id.o -cwd $index_command
    
done
    



