#!/bin/bash

set -u
set -e

if [ "$1" == -h ]; then
    echo ""
    echo "Maps reads from a fastq files to reference genome. Sorts and marks duplicate
with picard tools. Carries out indel realignment with GATK and runs various QC steps"
    echo ""
    echo "Paths to the  the reference genome file (REF_FILE),the path to the directory 
where mapped bam files (BAM_FILES) to be written and the path where the QC metrics folder 
(QC_METRICS) need to be specified inside the script"
    echo ""
    echo "The first command line argument is a 11 column tab-delimited text file with the
a paths to the two fastq files and followed by 9 columns with the read group infomation on 
each row with the format:"
    echo ""
    echo "fastq_1 fastq_2 rg_id rg_pl rg_pu rg_lb rg_pi rg_ds rg_dt rg_sm rg_cn"
    echo ""
    echo "example usage:"
    echo ""
    echo "sh bwa_mem_real.sh fastq_and_rg.txt"
    echo ""
    echo "Note: Line 71 needs altering if the fastq files are not gzipped"
    exit 1
fi


REF_FILE=/data/bo1pgc/gtit_1.04/genome/Parus_major_1.04.rename.fa # path to reference genome file that has been indexed 


#BGI files
                             
BAM_FILES=/fastdata/bo1pgc/parus_reseq/bgi_reseq/bams # path where bam files and intermediate file are to be written                
QC_METRICS=/data/bo1pgc/parus_reseq/bgi_reseq/mapping_qc #path where qc metrics will be written     
QSUB_PATH=/fastdata/bo1pgc/parus_reseq/bgi_reseq/qsub_out/bwa_mem_real # directory to write the qsub job log files
 


if [ ! -f "$REF_FILE.ann" ];  then
    echo "The reference fasta file needs to be indexed"
    exit 1
fi

if [ ! -d "$BAM_FILES" ]; then
    mkdir $BAM_FILES
fi

if [ ! -d "$QC_METRICS" ]; then
    mkdir $QC_METRICS
fi


if [ ! -d "$QSUB_PATH" ]; then
    mkdir $QSUB_PATH
fi


fastq_files=$1

cat $fastq_files | while read i j rg_id rg_pl rg_pu rg_lb rg_pi rg_ds rg_dt rg_sm rg_cn;

do
    bam_id="$(basename $i .clean.1.fq.gz)"
    
    fastq_1=$i
    fastq_2=$j

   
    TEMP=$QSUB_PATH/$bam_id 
    
    # create the temp directory needed by some of the picard tool steps
    if [ ! -d "$TEMP" ]; then
	mkdir $TEMP
    fi

    bwa_command="bwa mem -t 10 -M $REF_FILE $fastq_1 $fastq_2 | samtools view -b - > $BAM_FILES/$bam_id.bam"

    #bwa_command="bwa mem -t 20 -M $REF_FILE $fastq_1 $fastq_2 > $BAM_FILES/$bam_id.sam"
    
    # now this performs the sam to bam coversion by piping to samtools above
    #sam_bam="java -Xmx4g -jar ~/bin/picard.jar SamFormatConverter INPUT=$BAM_FILES/$bam_id.sam OUTPUT=$BAM_FILES/$bam_id.bam"
    #rm_sam="rm $BAM_FILES/$bam_id.sam"
    
    # Adds the read group info
    add_read_group="java -Xmx12g -jar ~/bin/picard.jar AddOrReplaceReadGroups I=$BAM_FILES/$bam_id.bam O=$BAM_FILES/$bam_id.id.bam \
RGID=$rg_id RGPL=$rg_pl RGPU=$rg_pu RGLB=$rg_lb RGPI=$rg_pi RGDS=$rg_ds RGDT=$rg_dt \
RGSM=$rg_sm RGCN=$rg_cn SORT_ORDER=coordinate TMP_DIR=$TEMP MAX_RECORDS_IN_RAM=2000000"

    sort_bam="java -Xmx12g -jar ~/bin/picard.jar SortSam I=$BAM_FILES/$bam_id.id.bam \
    O=$BAM_FILES/$bam_id.sorted.bam SORT_ORDER=coordinate TMP_DIR=$TEMP MAX_RECORDS_IN_RAM=2000000"
    
    #Mark duplicates
    mark_dups="java -Xmx64g -jar ~/bin/picard.jar MarkDuplicates INPUT=$BAM_FILES/$bam_id.sorted.bam \
  OUTPUT=$BAM_FILES/$bam_id.dedup.bam METRICS_FILE=$QC_METRICS/$bam_id.markdup.metrics.txt CREATE_INDEX=true TMP_DIR=$TEMP MAX_RECORDS_IN_RAM=12000000"
    
    #Perform indel realignment
    create_realign="java -Xmx8G -jar ~/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator \
 -I $BAM_FILES/$bam_id.dedup.bam -R $REF_FILE -o $BAM_FILES/$bam_id.dedup.bam.intervals"
    
    indel_realign="java -Xmx8G -jar ~/bin/GenomeAnalysisTK.jar -T IndelRealigner \
 -I $BAM_FILES/$bam_id.dedup.bam -R $REF_FILE -targetIntervals $BAM_FILES/$bam_id.dedup.bam.intervals \
 --consensusDeterminationModel USE_READS -o $BAM_FILES/$bam_id.dedup.real.bam"
 
    #Perform a QC step after the mark duplicate step 
    #run the validation step                                                                                                                        

    wgs_metrics="java -Xmx2g -jar ~/bin/picard.jar CollectWgsMetrics I=$BAM_FILES/$bam_id.dedup.real.bam\
 O=$QC_METRICS/$bam_id.wgsmetrics_file.txt R=$REF_FILE INCLUDE_BQ_HISTOGRAM=true"
    
    # run fastqc on mapped reads, reads that are not pcr duplicates and reads that are flagged primary alignments
    samtools_aligned="samtools view -b -F 1284 $BAM_FILES/$bam_id.dedup.real.bam > $BAM_FILES/$bam_id.aligned_reads.dedup.real.bam"
    fastqc_command="~/bin/FastQC/fastqc $BAM_FILES/$bam_id.aligned_reads.dedup.real.bam -o $QC_METRICS"
    
    validate_bam="java -Xmx4g -jar ~/bin/picard.jar ValidateSamFile I=$BAM_FILES/$bam_id.dedup.real.bam \
O=$QC_METRICS/$bam_id.validation.txt"

    #Cleanup
    
    rm_command="rm $BAM_FILES/$bam_id.id.bam $BAM_FILES/$bam_id.dedup.bai $BAM_FILES/$bam_id.aligned_reads.dedup.real.bam \
$BAM_FILES/$bam_id.sorted.bam $BAM_FILES/$bam_id.dedup.bam $BAM_FILES/$bam_id.dedup.bam.intervals $BAM_FILES/$bam_id.bam"
    
    
    #Submit the qsub job to icberg
    #Use -hold_jid to hold jobs until the required
    #previous job in the pipeline has been completed

    qsub -b y -q evolgen.q -P evolgen -pe openmp 10 -N bwa_mem.$bam_id -l mem=8G -l arch=intel* \
	-l h_rt=24:00:00 -e $QSUB_PATH/bwa_mem.$bam_id.e \
	-o $QSUB_PATH/bwa_mem.$bam_id.o -cwd $bwa_command
    
    #qsub -b y -q evolgen.q -P evolgen -N sam_to_bam.$bam_id -l mem=8G -l rmem=8G -l arch=intel* \
    #    -l h_rt=8:00:00 -hold_jid bwa_mem.$bam_id -e $QSUB_PATH/sam_to_bam.$bam_id.e \
    #	-o $QSUB_PATH/sam_to_bam.$bam_id.o -cwd $sam_bam
    
    #qsub -b y -q evolgen.q -P evolgen -N rm_sam.$bam_id -l mem=8G -l rmem=8G -l arch=intel* \
    #	-l h_rt=8:00:00 -hold_jid sam_to_bam.$bam_id -e $QSUB_PATH/rm_sam.$bam_id.e \
    #	-o $QSUB_PATH/rm_sam.$bam_id.o -cwd $rm_sam 

    qsub -b y -q evolgen.q -P evolgen -N add_rg.$bam_id -l mem=16G -l rmem=16G -l arch=intel* \
        -l h_rt=8:00:00 -hold_jid bwa_mem.* -e $QSUB_PATH/add_rg.$bam_id.e \
        -o $QSUB_PATH/add_rg.$bam_id.o -cwd $add_read_group
    
    qsub -b y -q evolgen.q -P evolgen -N sort_bams.$bam_id -l mem=16G -l rmem=16G -l arch=intel* \
        -l h_rt=8:00:00 -hold_jid add_rg.$bam_id -e $QSUB_PATH/sort_bams.$bam_id.e \
        -o $QSUB_PATH/sort_bams.$bam_id.o -cwd $sort_bam
    
    qsub -b y -q evolgen.q -P evolgen -N bam_dedup.$bam_id -l mem=100G -l rmem=100G  -l arch=intel* \
        -l h_rt=24:00:00 -hold_jid sort_bams.$bam_id -e $QSUB_PATH/bam_dedup.$bam_id.e \
        -o $QSUB_PATH/bam_dedup.$bam_id.o -cwd $mark_dups

    qsub -b y -q evolgen.q -P evolgen -N indel_target.$bam_id -l mem=12G -l rmem=12G  -l arch=intel* \
	-l h_rt=24:00:00 -hold_jid bam_dedup.$bam_id -e $QSUB_PATH/indel_target.$bam_id.e \
        -o $QSUB_PATH/indel_target.$bam_id.o -cwd $create_realign

    qsub -b y -q evolgen.q -P evolgen -N indel_realigner.$bam_id  -l mem=12G -l rmem=12G  -l arch=intel* \
	-l h_rt=24:00:00 -hold_jid indel_target.$bam_id -e $QSUB_PATH/indel_realigner.$bam_id.e \
        -o $QSUB_PATH/indel_realigner.$bam_id.o -cwd $indel_realign
    
#    qsub -b y -q evolgen.q -P evolgen -N validate_bam.$bam_id -l mem=10G -l rmem=10G  -l arch=intel* \
#	-l h_rt=8:00:00 -hold_jid indel_realigner.$bam_id -e $QSUB_PATH/validate_bam.$bam_id.e \
#	-o $QSUB_PATH/validate_bam.$bam_id.o -cwd $validate_bam 
    
#    qsub -b y -N wgs_metric.$bam_id -l mem=8G -l rmem=8G  -l arch=intel* \
#	-l h_rt=8:00:00 -hold_jid indel_realigner.$bam_id -e $QSUB_PATH/wgs_metric.$bam_id.e \
#	-o $QSUB_PATH/wgs_metric.$bam_id.o -cwd $wgs_metrics
    
#    qsub -b y -N samtools_aligned.$bam_id -l mem=8G -l rmem=8G  -l arch=intel* \
#	-l h_rt=8:00:00 -hold_jid indel_realigner.$bam_id -e $QSUB_PATH/samtools_aligned.$bam_id.e \
#	-o $QSUB_PATH/samtools_aligned.$bam_id.o -cwd $samtools_aligned
   
#    qsub -b y -N fastqc_bam.$bam_id -l mem=4G -l rmem=4G -l arch=intel* \
#        -l h_rt=8:00:00 -hold_jid samtools_aligned.$bam_id -e $QSUB_PATH/fastqc.$bam_id.e \
#        -o $QSUB_PATH/fastqc.$bam_id.o -cwd $fastqc_command
    
    qsub -b y -q evolgen.q -P evolgen -N remove_bams.$bam_id -l mem=8G -l rmem=8G -l arch=intel* \
        -l h_rt=8:00:00 -hold_jid indel_target.$bam_id -e $QSUB_PATH/remove_bams.$bam_id.e \
        -o $QSUB_PATH/remove_bams.$bam_id.o -cwd $rm_command
    
done
