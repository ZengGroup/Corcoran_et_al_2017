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
the genome on each line, reflecting the order of the
entries in the reference genome. This order is important 
as this sort order is assumed by the CatVariants tool run in the 
gathering stage. The scaffolds.intervals file must be in the directory 
where this script resides and must be named scaffolds.intervals"
    echo ""
    echo "Format of bam file list file:"
    echo ""
    echo "/path/to/bam<TAB>M"
    echo ""
    echo "Usage:"
    echo "sh haplotype_caller.sh bamfile_list.list chromosomes.txt"
    echo ""
    exit 1
fi

#Paths for calling with unrecalibrated bam files
REF_FILE=/data/bo1pgc/zebra_finch_genome/taeGut1.reordered.fa

#Paths for calling with recalibrated bam filess  
QSUB_PATH=/fastdata/bo1pgc/parus_reseq/singhal_zebra_finch_data/zf_snp_calling/qsub_out/haplotype_caller
GVCF_FILE=/fastdata/bo1pgc/parus_reseq/singhal_zebra_finch_data/zf_snp_calling/gVCF


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

bam_array=($(cat $1))
chromosome_array=($(cat $2)) 

for i in ${bam_array[@]}
do  
    sample_id="$(basename $i .recal.bam)" 
    
    for chrom in ${chromosome_array[@]} # scatter the calling across chromosomes
    do
		hap_caller="java -Xmx3g -jar ~/bin/GenomeAnalysisTK.jar -nct 2
						-T HaplotypeCaller
						-R $REF_FILE
						-I $i
						--emitRefConfidence GVCF
						-L $chrom
						-o $GVCF_FILE/$sample_id.$chrom.g.vcf"

		qsub -b y -q evolgen.q -P evolgen -pe openmp 2 -N hap_call.$sample_id.$chrom -l mem=6G -l rmem=6G  -l arch=intel* \
			-hold_jid print_reads.* -l h_rt=48:00:00 -e $QSUB_PATH/hap_call.$sample_id.$chrom.e \
			-o $QSUB_PATH/hap_call.$sample_id.$chrom.o -cwd $hap_caller
        	
    done

    #concatenate_script=/home/bo1pgc/parus_reseq/bgi_data/bgi_pipeline/concatenate_vcfs.py
# gather the chromosomal vcf files into one vcf for each sample
    python concatenate_vcfs.py -v gvcf -c $2 -l $QSUB_PATH -R $REF_FILE -p $GVCF_FILE -q general -n $sample_id -o $sample_id.raw.snps.indels.g.vcf.gz

done
