#!/bin/bash

#$ -l h_rt=08:00:00                                                                                                                                                              
#$ -l mem=8G                                                                                                                                                                     
#$ -l rmem=8G                                                                                                                                                                    
#$ -l arch=intel*                                                                                                                                                                
#$ -e analyse_covariates.e                                                                                                                                                       
#$ -o analyse_covariates.o



if [ "$1" == -h ]; then
    echo ""
    echo "Runs the BQSR analyse covariates step from the GATK best practice pipeline
to produce the before and after BQSR plots.

The first command line argument is a text file listing the realigned
full paths to unrecalibrated BAM files in the first column and the 
sex of the birds in the second column. The sample names are parse from this text file"
    echo ""
    echo "Format of file:"
    echo ""
    echo "/path/to/bam<TAB>M"
    echo ""
    echo "Usage:"
    echo "sh analyse_covariate.sh bamfile_list.txt  "
    echo ""
    exit 1
fi

# PATHS  required files and output directories for files produced


QSUB_PATH=/fastdata/bo1pgc/parus_reseq/bgi_reseq/qsub_out/bqsr
REF_FILE=/data/bo1pgc/gtit_1.04/genome/Parus_major_1.04.rename.fa
BQSR_REPORTS=/fastdata/bo1pgc/parus_reseq/bgi_reseq/bqsr_reports


if [ ! -f "$REF_FILE" ]; then
    echo "The path to the reference file is incorrect"
    exit 1
fi

if [ ! -d "SBQSR_REPORTS" ]; then
    mkdir $BQSR_REPORTS
fi


if [ ! -d "$QSUB_PATH" ]; then
    mkdir $QSUB_PATH
fi


module load apps/R/3.2.1

cut -f1 $1 | while read i ;
do 
    bam_id="$(basename $i .dedup.real.bam)"

    java -Xmx4g -jar ~/bin/GenomeAnalysisTK.jar -T AnalyzeCovariates  -R $REF_FILE \
	-before $BQSR_REPORTS/$bam_id.recal_data.table \
	-after $BQSR_REPORTS/$bam_id.post_recal_data.table \
	-plots $BQSR_REPORTS/$bam_id.recalibration_plots.pdf    
done