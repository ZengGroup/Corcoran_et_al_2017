#!/bin/bash                                            

set -e
set -u                                                                       

if [ "$1" == -h ]; then
    echo ""
    echo "Applies the GATK recommended hard filters to a vcf
containing SNPs and INDELs. This script creates a two separate 
SNP and INDEL filtered vcf files. The hard filters appied to 
vcf files are taken from the GATK documentation at 
http://gatkforums.broadinstitute.org/discussion/2806/howto-apply-hard-filters-to-a-call-set"
    echo ""
    echo ""
    echo "Usage:"
    echo "sh hard_filter_vcf.sh vcf_file"
    echo ""
    exit 1
fi

#PATHS
QSUB_PATH=/fastdata/bo1pgc/parus_reseq/bgi_reseq/qsub_out/hard_filter_vcf
REF_FILE=/data/bo1pgc/gtit_1.04/genome/Parus_major_1.04.rename.fa
VCF_FILE=/fastdata/bo1pgc/parus_reseq/bgi_reseq/vcfs_for_bqsr # directory to ouput filtered vcfs

#queue="-q evolgen.q -P evolgen"
queue=""

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


output_name="$(basename $1 .raw.snps.indels.vcf.gz)"

# Select SNPs
select_snps="java -Xmx8G -jar ~/bin/GenomeAnalysisTK.jar -T SelectVariants  
    -R $REF_FILE -V $1 
    -selectType SNP -o $VCF_FILE/$output_name.raw_snps.vcf"

#Select INDELs
select_indels="java -Xmx8G -jar ~/bin/GenomeAnalysisTK.jar -T SelectVariants 
    -R $REF_FILE -V $1
    -selectType INDEL -o $VCF_FILE/$output_name.raw_indels.vcf" 

#Apply GATK Hard FIlters
filter_snps="java -Xmx8G -jar ~/bin/GenomeAnalysisTK.jar -T VariantFiltration  
    -R $REF_FILE  
    -V $VCF_FILE/$output_name.raw_snps.vcf 
    --filterExpression \"(vc.hasAttribute('QD') && QD < 2.0) || FS > 60.0 || SOR > 4.0 || MQ < 40.0 || 
(vc.hasAttribute('MQRankSum') && MQRankSum < -12.5) || (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0)\"  
    --filterName \"gatk_hard_snp\" 
    -o $VCF_FILE/$output_name.filtered_snps.vcf"

filter_indels="java -Xmx8G -jar ~/bin/GenomeAnalysisTK.jar -T VariantFiltration  
    -R $REF_FILE  
    -V $VCF_FILE/$output_name.raw_indels.vcf
    --filterExpression \"(vc.hasAttribute('QD') && QD < 2.0) || FS > 200.0 || (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0)\" 
    --filterName \"gatk_hard_indel\"  
    -o $VCF_FILE/$output_name.filtered_indels.vcf"

# create filtered vcf files
exclude_filter_snps="java -Xmx8G -jar ~/bin/GenomeAnalysisTK.jar -T SelectVariants      
    -R $REF_FILE -V $VCF_FILE/$output_name.filtered_snps.vcf                                                       
    -o $VCF_FILE/$output_name.gatk_hard_filter.snp.vcf.gz --excludeFiltered"

exclude_filter_indels="java -Xmx8G -jar ~/bin/GenomeAnalysisTK.jar -T SelectVariants
    -R $REF_FILE -V $VCF_FILE/$output_name.filtered_indels.vcf                                                                                                              
    -o $VCF_FILE/$output_name.gatk_hard_filter.indel.vcf.gz --excludeFiltered"

# remove intermediate vcf files that were generated
rm_temp="rm $VCF_FILE/$output_name.filtered_indels.vcf $VCF_FILE/$output_name.filtered_snps.vcf $VCF_FILE/$output_name.raw_indels.vcf
$VCF_FILE/$output_name.raw_snps.vcf"


qsub -b y $queue -N select_snps.$output_name -l mem=12G -l rmem=12G  -l arch=intel* \
          -l h_rt=8:00:00 -e $QSUB_PATH/select_snps.$output_name.e \
          -o $QSUB_PATH/select_snps.$output_name.o -cwd $select_snps

qsub -b y $queue -N select_indels.$output_name -l mem=12G -l rmem=12G  -l arch=intel* \
          -l h_rt=8:00:00 -e $QSUB_PATH/filter_indels.$output_name.e \
          -o $QSUB_PATH/select_indels.$output_name.o -cwd $select_indels

qsub -b y $queue -N filter_snps.$output_name -l mem=12G -l rmem=12G  -l arch=intel* \
          -l h_rt=8:00:00 -hold_jid select_snps.$output_name -e $QSUB_PATH/filter_snps.$output_name.e \
          -o $QSUB_PATH/filter_snps.$output_name.o -cwd $filter_snps

qsub -b y $queue -N filter_indels.$output_name -l mem=12G -l rmem=12G  -l arch=intel* \
          -l h_rt=8:00:00 -hold_jid select_indels.$output_name -e $QSUB_PATH/filter_indels.$output_name.e \
          -o $QSUB_PATH/filter_indels.$output_name.o -cwd $filter_indels

qsub -b y $queue -N exclude_filter_snps.$output_name -l mem=12G -l rmem=12G  -l arch=intel* \
          -l h_rt=8:00:00 -hold_jid filter_snps.$output_name -e $QSUB_PATH/exclude_filter_snps.$output_name.e \
          -o $QSUB_PATH/exclude_filter_snps.$output_name.o -cwd $exclude_filter_snps

qsub -b y $queue -N exclude_filter_indels.$output_name -l mem=12G -l rmem=12G  -l arch=intel* \
          -l h_rt=8:00:00 -hold_jid filter_indels.$output_name -e $QSUB_PATH/exclude_filter_indels.$output_name.e \
          -o $QSUB_PATH/exclude_filter_indels.$output_name.o -cwd $exclude_filter_indels

qsub -b y $queue -N rm_temp.$output_name -l mem=12G -l rmem=12G  -l arch=intel* \
          -l h_rt=8:00:00 -hold_jid exclude_filter_*.* -e $QSUB_PATH/rm_temp.$output_name.e \
          -o $QSUB_PATH/rm_temp.$output_name.o -cwd $rm_temp