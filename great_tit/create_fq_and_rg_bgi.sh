#!/bin/bash

FASTQ_PATH=/fastdata/bo1pgc/parus_reseq/bgi_reseq/fastq_files # path to fastq files
FASTQ_CLEAN_PATH=/fastdata/bo1pgc/parus_reseq/bgi_reseq/clean_fastq_files #pat  where clean files will be written

ls $FASTQ_PATH/*/clean_reads/*.fq.gz \
 | cut -f 7,9 -d'/' \
 | python get_read_group_info_bgi.py $FASTQ_PATH $FASTQ_CLEAN_PATH
