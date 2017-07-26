#!/bin/bash

set -e
set -u


if [ "$1" == -h ]; then
    echo ""
    echo "Runs the GenotypeGVCFs tool on a listof gvcf files, produced by 
the haplotype caller,in a scatter and gather approach. 
The first command line argument is a text file listing 
the path to the samples' gVCF files to run. 

The second command line argument lists the chromosome names in 
the genome on each line, reflecting the order of the
entries in the reference genome fasta file. This order is important 
as this sort order is assumed by the CatVariants tool run in the 
gathering stage. The scaffolds.intervals file must be in the directory 
where this script resides and must be named scaffolds.intervals.

The third command line argument is the name that will be prepended to the
output vcf file that will be produced. The ouput file produce will have the
name outfile_name.raw.snps.indels.vcf

Paths to the sourve vcf files, the reference genome file and the path where
log file should be written need to be specified within the script."
    echo ""
    echo ""
    echo "Usage:"
    echo "sh genotype_gvcfs.sh gvcf_list.list chromosomes.txt outfile_name"
    echo ""
    exit 1
fi


# analysis with base recal gvcf
REF_FILE=/data/bo1pgc/zebra_finch_genome/taeGut1.reordered.fa

#Paths for calling with recalibrated bam filess  
QSUB_PATH=/fastdata/bo1pgc/parus_reseq/singhal_zebra_finch_data/zf_snp_calling/qsub_out/genotype_gcvf
VCF_FILE=/fastdata/bo1pgc/parus_reseq/singhal_zebra_finch_data/zf_snp_calling/VCF

queue="-q evolgen.q -P evolgen"

if [ ! -f "$REF_FILE" ]; then
    echo "The path to the reference file is incorrect"
    exit 1
fi

if [ ! -d "$QSUB_PATH" ]; then
    mkdir $QSUB_PATH
fi

if [ ! -d "$VCF_FILE" ]; then
    mkdir $VCF_FILE
fi


gvcf_list=$1
chromosome_array=($(cat $2)) 
outfile_name=$3


for chrom in ${chromosome_array[@]} # scatter the calling across chromosomes
do

    genotyper="java -Xmx6g -jar ~/bin/GenomeAnalysisTK.jar 
               --disable_auto_index_creation_and_locking_when_reading_rods                                                  
               -T GenotypeGVCFs -allSites                                                                                                                
               -R $REF_FILE                                                                                                                                          
               -L $chrom
               -V $gvcf_list
               -hets 0.006                                                                                                                      
               -o $VCF_FILE/$outfile_name.$chrom.raw.snps.indels.vcf"
	   
    qsub -b y $queue -N genotype.$outfile_name.$chrom -l mem=10G -l rmem=10G -l arch=intel* \
    	 -l h_rt=24:00:00 -e $QSUB_PATH/genotype.$outfile_name.$chrom.e \
     	 -o $QSUB_PATH/genotype.$outfile_name.$chrom.o -cwd $genotyper
       
done

# gather the chromosomal vcf files into one vcf for the whole genome 
python concatenate_vcfs.py -v vcf -c $2 -l $QSUB_PATH -R $REF_FILE -p $VCF_FILE -n $outfile_name -o $outfile_name.raw.snps.vcf.all_sites.temp.vcf.gz -q general

remove_alternate="java -Xmx3g -jar ~/bin/GenomeAnalysisTK.jar 
                 -T SelectVariants 
                 -R $REF_FILE 
                 -V $VCF_FILE/$outfile_name.raw.snps.vcf.all_sites.temp.vcf.gz 
                 -trimAlternates 
                 -o $VCF_FILE/$outfile_name.raw.snps.all_sites.vcf.gz"

qsub -b y $queue -N remove_alt.$outfile_name -l mem=10G -l rmem=10G -l arch=intel* \
     -hold_jid gather_vcf.*  -l h_rt=24:00:00 -e $QSUB_PATH/remove_alternate.$outfile_name.e \
     -o $QSUB_PATH/remove_alternate.$outfile_name.o -cwd $remove_alternate