#!/bin/bash

set -e
set -u

#$ -l h_rt=8:00:00
#$ -l mem=40G
#$ -l rmem=40G
#$ -l arch=intel*
#$ -e /fastdata/bo1tg/err/zfinch_vqsr.e
#$ -o /fastdata/bo1tg/out/zfinch_vqsr.o
#$ -N zf_vqsr
#$ -q evolgen.q -P evolgen

module load apps/binapps/GATK
module load apps/R/3.2.1 # need to install ggplot2

# PATHS required files and output directories for files produced
REF_FILE=/data/bo1tg/zebra/taeGut1.reordered.fa
VQSR_REPORTS=/fastdata/bo1tg/VCF_10birds/VCF_merge/vqsr_reports
TYPE=SNP #specify SNP or INDEL
VCF_FILE=/fastdata/bo1tg/VCF_10birds/zf_10birds.raw.snps.all_sites.vcf.gz
TRAINING_VCF=/fastdata/bo1tg/VCF_10birds/VCF_merge/Merged_VCF/0002sorted.vcf # includes Z chromosome
NAME=gt_10birds

if [ "$TYPE" != "SNP" -a "$TYPE" != "INDEL" ]; then
echo "Need to specify SNP or INDEL as first argument"
exit 1
fi

if [ ! -f "$REF_FILE" ]; then
echo "The path to the reference file is incorrect"
exit 1
fi

if [ ! -d "$VQSR_REPORTS" ]; then
mkdir $VQSR_REPORTS
fi

if [ ! -f "$VCF_FILE" ]; then
echo "The path to the raw vcf file is incorrect"
exit 1
fi

if [ ! -f "$TRAINING_VCF" ]; then
echo "The path to the training VCF file is incorrect"
exit 1
fi

java -Xmx24g -jar $GATKHOME/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_FILE \
-input $VCF_FILE \
-resource:zf_reseq,known=true,training=true,truth=true,prior=10.0 $TRAINING_VCF \
-an DP \
-an QD \
-an FS \
-an SOR \
-an MQRankSum \
-an ReadPosRankSum \
-an InbreedingCoeff \
-mode $TYPE \
-tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 90.0 \
-recalFile $VQSR_REPORTS/$NAME.recalibrate_SNP.recal \
-tranchesFile $VQSR_REPORTS/$NAME.recalibrate_SNP.tranches \
-rscriptFile $VQSR_REPORTS/$NAME.recalibrate_SNP_plots.R