#!/bin/bash
set -e
set -u
#$ -l h_rt=24:00:00
#$ -l mem=4G
#$ -l rmem=4G
#$ -l arch=intel*
#$ -e /fastdata/bo1tg/err/sort_bgzip.e
#$ -o /fastdata/bo1tg/out/sort_bgzip.o
#$ -N sort_bgzip
#$ -q evolgen.q -P evolgen


infileA=/fastdata/bo1tg/VCF_10birds/VCF_merge/gt_10birds.gatk.SNPsONLY.vcf
outfileA=/fastdata/bo1tg/VCF_10birds/VCF_merge/gt_10birds.gatk.SNPsONLY.sorted.vcf
cat $infileA | vcf-sort -t /fastdata/bo1tg/tmp > $outfileA
bgzip $outfileA 

infileB=/fastdata/bo1tg/VCF_10birds/VCF_merge/gt_10birds.freebayes.SNPsONLY.vcf
outfileB=/fastdata/bo1tg/VCF_10birds/VCF_merge/gt_10birds.freebayes.SNPsONLY.sorted.vcf
cat $infileB | vcf-sort -t /fastdata/bo1tg/tmp > $outfileB
bgzip $outfileB

tabix -p vcf $outfileB.gz
tabix -p vcf $outfileA.gz
