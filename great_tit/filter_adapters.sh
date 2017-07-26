#!/bin/bash

set -e
set -u

if [ "$1" == -h ]; then
    echo ""
    echo "Remove adapter contamination with the program CutAdapt 
and run fastqc on the cleaned fastq files"
    echo ""
    echo "Usage:"   
    echo ""
    echo "sh filter_adapters.sh fastq_and_rg_info.txt"
    exit 1
fi



QSUB_PATH=/fastdata/bo1pgc/parus_reseq/bgi_reseq/qsub_out/adapter_removal
CLEAN_FASTQ=/fastdata/bo1pgc/parus_reseq/bgi_reseq/clean_fastq_files
OUTDIR=/data/bo1pgc/parus_reseq/bgi_reseq/fastqc_cutadapt_cleaned


if [ ! -d "$QSUB_PATH" ]; then
    mkdir $QSUB_PATH
fi

if [ ! -d "$CLEAN_FASTQ" ]; then
    mkdir $CLEAN_FASTQ
fi

if [ ! -d "$OUTDIR" ]; then
    mkdir $OUTDIR
fi



# adapters for sample 917 BGI data 
adapter_1_917="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"
adapter_2_917="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNNNNNAGATCTCGGTGGTCGCCGTA"
five_prime_1_917="CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT" # rev_comp adapter 1
five_prime_2_917="AATGATACGGCGACCACCGAGATCTNNNNNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT" # rev_comp adapter 2
five_prime_seq1_917="CGTCATA" # overrespresented in fq1 of 917
five_prime_seq2_917="CGGTGGT" # this sequence in found on 5' from position 1 at  ~29000 reads in fq2 of 917 sample


#adapters sequeunce to trim for the 9 other BGI samples
adapter_1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACACGCGGATCTCGTATGCCGTCTTCTGCTTG"
adapter_2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
five_prime_1="CAAGCAGAAGACGGCATACGAGATCCGCGTGTGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
five_prime_2="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"

five_prime_seq1="CTACACGACGCTCTTCCGATCT"  # enriched from pos 12-29 based on fastqc kmer analysis                                      
five_prime_seq2="TACACTCTT" # enriched in fg1 and fq2 from pos 1 analysis

five_prime_seq3="CGTCATA" #enriched from pos 1 in fq1
five_prime_seq4="NNNNNNGCGTCGT" # enriched from pos 7 in fq1
five_prime_seq5="NNNCTAGACG" # enriched from pos 4 in fq1 and fq2
five_prime_seq6="TCAGACG" # enriced from pos 1 in fq2
five_prime_seq7="TCAGACGTGTGCTCTTCCGATC" # enriched from pos 1 in TR43666 sample fq2 afer frist round of trimming
five_prime_seq8="GATCGGA" # enriched from pos 1 of fq2 in TR43666, 15 and 167 



cut -f1,2,10 $1 | while read fq1 fq2 sample; # loop for bgi data

do
    
    clean_file=$sample # rename bgi fastq files to have sample name
    
    if [ $sample == '917' ]; then # bgi 917 filtering
	
	filter_adapters="~/.local/bin/cutadapt -m 70 -n 2 -b $adapter_1_917 -B $adapter_2_917 -g $five_prime_1_917 -g ^$five_prime_seq1_917 -G $five_prime_2_917 -G ^$five_prime_seq2_917 -o $CLEAN_FASTQ/$clean_file.clean.1.fq.gz -p $CLEAN_FASTQ/$clean_file.clean.2.fq.gz $fq1 $fq2"
    
	
    else # bgi 9 samples filtering

	filter_adapters="~/.local/bin/cutadapt -m 70 -n 2 -b $adapter_1 -B $adapter_2 -g $five_prime_1 -g $five_prime_2 -g $five_prime_seq1 -g ^$five_prime_seq2 -g ^$five_prime_seq3 -g ^$five_prime_seq4 -g ^$five_prime_seq5  -G $five_prime_2 -G $five_prime_seq1 -G ^$five_prime_seq2 -G ^$five_prime_seq5 -G ^$five_prime_seq6 -G ^$five_prime_seq7 -G ^$five_prime_seq8 -o $CLEAN_FASTQ/$clean_file.clean.1.fq.gz -p $CLEAN_FASTQ/$clean_file.clean.2.fq.gz $fq1 $fq2"

    fi
    
    echo $filter_adapters

    fastqc="~/bin/FastQC/fastqc -t 2 $CLEAN_FASTQ/$clean_file.clean.1.fq.gz $CLEAN_FASTQ/$clean_file.clean.2.fq.gz -o $OUTDIR"
    
    queue="-q evolgen.q -P evolgen"
    #queue=""
    
    qsub -b y $queue -N filter_adapter.$clean_file -l mem=8G -l rmem=8G -l arch=intel* \
	-l h_rt=48:00:00 -e $QSUB_PATH/filter_adapters.$clean_file.e \
        -o $QSUB_PATH/filter_adapters.$clean_file.o -cwd $filter_adapters
    
    qsub -b y $queue -pe openmp 2 -N fastqc_clean.$clean_file -l mem=8G -l rmem=8G -l arch=intel* \
        -l h_rt=8:00:00 -hold_jid filter_adapter.$clean_file -e $QSUB_PATH/fastqc_clean.$clean_file.e \
        -o $QSUB_PATH/fastqc_clean.$clean_file.o -cwd $fastqc

done