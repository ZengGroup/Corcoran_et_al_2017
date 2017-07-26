# run before module load apps/binapps/GATK
import os, sys

# extract tranches (take two 99.9 and 99 and 90 to test)

gatk_exe='java -Xmx8g -jar $GATKHOME/GenomeAnalysisTK.jar'
ref_fasta='/fastdata/bo1tg/ref_genome/Parus_major_1.04.rename.fa'
vcf='/fastdata/bo1tg/VCF_10birds/bgi_10birds.raw.snps.indels.all_sites.vcf.gz'
recal_file='/fastdata/bo1tg/VCF_10birds/VCF_merge/vqsr_reports/gt_10birds.recalibrate_SNP.recal'
tranch_file='/fastdata/bo1tg/VCF_10birds/VCF_merge/vqsr_reports/gt_10birds.recalibrate_SNP.tranches'


for tranch in [90,99,99.9,100]:
	output='/fastdata/bo1tg/VCF_10birds/VCF_merge/vqsr/recalibrated_snps_'+str(tranch)+'.vcf.gz'

	command=gatk_exe+' -T ApplyRecalibration -R '+ref_fasta+' -input '+vcf+' -mode SNP --ts_filter_level '+str(tranch)+' -recalFile '+recal_file+' -tranchesFile '+tranch_file+' -o '+output
	#print(gatk_exe+' -T ApplyRecalibration -R '+ref_fasta+' -input '+vcf+' -mode SNP --ts_filter_level '+str(tranch)+' -recalFile '+recal_file+' -tranchesFile '+tranch_file+' -o '+output)
	tranch=str(tranch)
	#os.system('qsub -q long.q -l rmem=20G -l mem=20G -b y -N rec'+tranch+' -e err/wSF'+str(tranch)+'.e -o out/wSF'+str(tranch)+'.o -cwd '+command)
	os.system('qsub -q evolgen.q -P evolgen -l rmem=20G -l mem=20G -b y -N rec'+tranch+' -e /fastdata/bo1tg/err/gt_ar'+str(tranch)+'.e -o /fastdata/bo1tg/out/gt_ar'+str(tranch)+'.o -cwd '+command)


