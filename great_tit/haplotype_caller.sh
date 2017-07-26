#!/bin/bash

set -e
set -u

if [ "$1" == -h ]; then
    echo ""
    echo "Runs the haplotype caller in GVCF mode on a list
of bam files in a scatter and gather approach. The first 
command line argument is a tab-delimited text file listing 
the path to each sample BAM file to run through the haplotype 
caller in column 1,  and the sex of the bird in column 2 (M or F). 

The second command line argument lists the chromosome names in 
the genome on each line and the scaffolds.intervalfile is listed 
as the second last entry before chrZ, reflecting the order of the
entries in the reference genome fasta file. This order is important 
as this sort order is assumed by the CatVariants tool run in the 
gathering stage. The scaffolds.intervals file must be in the directory 
where this script resides and must be named scaffolds.intervals"
    echo ""
    echo "Format of bam file list file:"
    echo ""
    echo "/path/to/bam<TAB>M"
    echo ""
    echo "Usage:"
    echo "sh haplotype_caller.sh bamfile_list.txt chromosomes.txt"
    echo ""
    exit 1
fi

#Paths for calling with unrecalibrated bam files
REF_FILE=/data/bo1pgc/gtit_1.04/genome/Parus_major_1.04.rename.fa


#QSUB_PATH=/fastdata/bo1pgc/parus_reseq/bgi_reseq/qsub_out/unrecal_gvcf
#GVCF_FILE=/fastdata/bo1pgc/parus_reseq/bgi_reseq/unrecal_gvcf

#Paths for calling with recalibrated bam filess  
QSUB_PATH=/fastdata/bo1pgc/parus_reseq/bgi_reseq/qsub_out/recal_gvcf                                                                                                         
GVCF_FILE=/fastdata/bo1pgc/parus_reseq/bgi_reseq/recal_gvcf

queue="-q evolgen.q -P evolgen" 

if [ ! -f "$REF_FILE" ]; then
    echo "The path to the reference file is incorrect"
    exit 1
fi

if [ ! -d "$QSUB_PATH" ]; then
    mkdir $QSUB_PATH
fi

if [ ! -d "$GVCF_FILE" ]; then
    mkdir $GVCF_FILE
fi


chromosome_array=($(cat $2)) 

num=0
cat $1 | while read i j
do  
    sample_id="$(basename $i .dedup.real.*)" 
    
    for chrom in ${chromosome_array[@]} # scatter the calling across chromosomes
    do
	
	if [ $chrom == 'chrZ' -a $j == 'F' ]; then # check for females
	    
	    hap_caller="java -Xmx3g -jar ~/bin/GenomeAnalysisTK.jar -nct 4
                       -T HaplotypeCaller 
                       -R $REF_FILE 
                       -I $i 
                       --emitRefConfidence GVCF 
                       -L $chrom 
                       -ploidy 1 
                       -o $GVCF_FILE/$sample_id.$chrom.raw.snps.indels.g.vcf"
	
	    
	else

	    if [ $chrom == 'scaffolds.intervals' ]; then
               chrom='scaffolds'
	       
	       hap_caller="java -Xmx3g -jar ~/bin/GenomeAnalysisTK.jar -nct 4
                           -T HaplotypeCaller 
                           -R $REF_FILE
                           -I $i
                           --emitRefConfidence GVCF
                           -L $chrom.intervals
                           -o $GVCF_FILE/$sample_id.$chrom.raw.snps.indels.g.vcf"

	    else

		hap_caller="java -Xmx3g -jar ~/bin/GenomeAnalysisTK.jar -nct 4
                           -T HaplotypeCaller                                                                                
                           -R $REF_FILE                                                                                      
                           -I $i                                                                                             
                           --emitRefConfidence GVCF                                                                          
                           -L $chrom                                                                      
                           -o $GVCF_FILE/$sample_id.$chrom.raw.snps.indels.g.vcf"
	
	    fi
	
	fi
	
	
	qsub -b y $queue -pe openmp 4 -N hap_call.$sample_id.$chrom -l mem=6G -l rmem=6G  -l arch=intel* \
	    -hold_jid print_reads.* -l h_rt=24:00:00 -e $QSUB_PATH/hap_call.$sample_id.$chrom.e \
	    -o $QSUB_PATH/hap_call.$sample_id.$chrom.o -cwd $hap_caller
        
	
    done

# gather the chromosomal vcf files into one vcf for each sample
    python concatenate_vcfs.py -v gvcf -c $2 -l $QSUB_PATH -R $REF_FILE -p $GVCF_FILE -q general -n $sample_id  -o $sample_id.raw.snps.indels.g.vcf.gz


done
