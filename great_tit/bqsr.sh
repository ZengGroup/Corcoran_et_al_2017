#!/bin/bash

if [ "$1" == -h ]; then
    echo ""
    echo "Runs the BQSR from the GATK best practice pipeline
on multiple unrecalibrated BAM file listed in a file. The list 
of bam files is given as commandline argument 1. The hard filtered VCF files 
produced from the unrecalibrated bam files are used as the known sites 
in this analysis. The paths to reference genome si given within the script 
and the VCF files needed are given at commandline arguments 2 and 3 . 
The path to where the BQSR reports and the recalibrated BAM files 
are to be written must also be specified within the script.

The first command line argument is a text file listing the realigned
unrecalibrated BAM files that will be recalibrated."
    echo ""
    echo "Format of file:"
    echo ""
    echo "/path/to/bam<TAB>M"
    echo ""
    echo "Usage:"
    echo "sh bqsr.sh bamfile_list.txt snp.vcf indel.vcf "
    echo ""
    exit 1
fi

# PATHS  required files and output directories for files produced

BAM_FILES=/fastdata/bo1pgc/parus_reseq/bgi_reseq/recalibrated_bams
QSUB_PATH=/fastdata/bo1pgc/parus_reseq/bgi_reseq/qsub_out/bqsr
REF_FILE=/data/bo1pgc/gtit_1.04/genome/Parus_major_1.04.rename.fa
BQSR_REPORTS=/fastdata/bo1pgc/parus_reseq/bgi_reseq/bqsr_reports

#queue="-q evolgen.q -P evolgen"
queue=""

if [ ! -d "$BAM_FILES" ]; then
    mkdir $BAM_FILES
fi

if [ ! -f "$REF_FILE" ]; then
    echo "The path to the reference file is incorrect"
    exit 1
fi

if [ ! -d "SBQSR_REPORTS" ]; then
    mkdir $BQSR_REPORTS
fi

if [ ! -f "$2" ]; then
    echo "The path to the SNP vcf file is incorrect"
    exit 1
fi

if [ ! -f "$3" ]; then
    echo "The path to the INDEL VCF file is incorrect"
    exit 1
fi

if [ ! -d "$QSUB_PATH" ]; then
    mkdir $QSUB_PATH
fi

VCF_SNP=$2
VCF_INDEL=$3

TEMP=$QSUB_PATH/temp

if [ ! -d "$TEMP" ]; then
    mkdir $TEMP
fi



cut -f1 $1 | while read i ;
do 
    bam_id="$(basename $i .dedup.real.bam)"
    
#Analyze patterns of covariation in the sequence dataset
    bqsr_1="java -Xmx9g -jar ~/bin/GenomeAnalysisTK.jar -T BaseRecalibrator --disable_auto_index_creation_and_locking_when_reading_rods 
    -R $REF_FILE -nct 2 
    -I $i  
    -knownSites $VCF_SNP 
    -knownSites $VCF_INDEL 
    -o $BQSR_REPORTS/$bam_id.recal_data.table" 

##analyze covariation remaining after recalibration
#    bqsr_2="java -Xmx9g -jar ~/bin/GenomeAnalysisTK.jar -T BaseRecalibrator --disable_auto_index_creation_and_locking_when_reading_rods 
#    -R $REF_FILE -nct 2
#    -I $i  
#    -knownSites $VCF_SNP 
#    -knownSites $VCF_INDEL  
#    -BQSR $BQSR_REPORTS/$bam_id.recal_data.table 
#    -o $BQSR_REPORTS/$bam_id.post_recal_data.table "

Generate before/after plots
    analyse_covariates="module load apps/R/3.2.1; java -Xmx8g -jar ~/bin/GenomeAnalysisTK.jar -T AnalyzeCovariates 
    -R $REF_FILE
    -before $BQSR_REPORTS/$bam_id.recal_data.table
    -after $BQSR_REPORTS/$bam_id.post_recal_data.table 
    -plots $BQSR_REPORTS/$bam_id.recalibration_plots.pdf"

# apply recalibration to reads
    print_reads="java -Djava.io.tmpdir=$TEMP -Xmx18g -jar ~/bin/GenomeAnalysisTK.jar -T PrintReads 
    -R $REF_FILE
    -I $i 
    -BQSR $BQSR_REPORTS/$bam_id.recal_data.table 
    -o $BAM_FILES/$bam_id.dedup.real.recal.bam"

    qsub -b y -pe openmp 2 $queue -N bqsr_1.$bam_id -l mem=12G -l rmem=12G  -l arch=intel* \
	-l h_rt=48:00:00 -hold_jid exclude_filter_*.* -e $QSUB_PATH/bqsr_1.$bam_id.e \
	-o $QSUB_PATH/bqsr_1.$bam_id.o -cwd $bqsr_1

    qsub -b y -pe openmp 2 $queue -N bqsr_2.$bam_id -l mem=12G -l rmem=12G  -l arch=intel* \
        -l h_rt=48:00:00 -hold_jid bqsr_1.$bam_id  -e $QSUB_PATH/bqsr_2.$bam_id.e \
        -o $QSUB_PATH/bqsr_2.$bam_id.o -cwd $bqsr_2
    
#    qsub -b y $queue -N analyse_covariates.$bam_id -l mem=12G -l rmem=12G  -l arch=intel* \
#        -l h_rt=48:00:00 -hold_jid bqsr_2.$bam_id -e $QSUB_PATH/analyse_covariates.$bam_id.e \
#        -o $QSUB_PATH/analyse_covariates.$bam_id.o -cwd $analyse_covariates
   
   
    qsub -b y -N print_reads.$bam_id -l mem=24G -l rmem=24G -l arch=intel* \
        -l h_rt=48:00:00 -hold_jid bqsr_1.$bam_id -e $QSUB_PATH/print_reads.$bam_id.e \
        -o $QSUB_PATH/print_reads.$bam_id.o -cwd $print_reads

done