import os

infile_gatk='/fastdata/bo1tg/VCF_10birds/bgi_10birds.raw.snps.indels.all_sites.vcf.gz'
infile_freebayes='/fastdata/bo1tg/VCF_10birds/VCF_merge/freebayes_merged.vcf'


# extract SNPs only

minD=25.
maxD=100.
minQ=20.

command='vcftools --gzvcf '+infile_gatk+' --remove-indels --recode --recode-INFO-all --max-missing 1 --maf 0.000001 --out SNPs_only --minQ '+str(minQ)+' --min-meanDP '+str(minD)+' --max-meanDP '+str(maxD)+' -c'
os.system('qsub -q long.q -b y -N extract -e /fastdata/bo1tg/err/extract_vcf.e -o /fastdata/bo1tg/VCF_10birds/VCF_merge/gt_10birds.gatk.SNPsONLY.vcf -cwd '+command)
command='vcftools --vcf '+infile_freebayes+' --remove-indels --recode --recode-INFO-all --max-missing 1 --maf 0.000001 --out SNPs_only --minQ '+str(minQ)+' --min-meanDP '+str(minD)+' --max-meanDP '+str(maxD)+' -c'
os.system('qsub -q long.q -b y -N extract -e /fastdata/bo1tg/err/extract_vcf.e -o /fastdata/bo1tg/VCF_10birds/VCF_merge/gt_10birds.freebayes.SNPsONLY.vcf -cwd '+command)
