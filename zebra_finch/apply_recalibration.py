# run before module load apps/binapps/GATK
import os

gatk_exe = 'java -Xmx8g -jar $GATKHOME/GenomeAnalysisTK.jar'
ref_fasta = '/data/bo1tg/zebra/taeGut1.reordered.fa'
vcf = '/fastdata/bo1tg/zebra_VCF_10birds/zf_10birds.raw.snps.all_sites.vcf.gz'
recal_file = '/fastdata/bo1tg/zebra_VCF_10birds/VCF_merge/vqsr_reports/zf_10birds.recalibrate_SNP.recal'
tranch_file = '/fastdata/bo1tg/zebra_VCF_10birds/VCF_merge/vqsr_reports/zf_10birds.recalibrate_SNP.tranches'


for tranch in [90, 99, 99.9, 100]:
	output = '/fastdata/bo1tg/zebra_VCF_10birds/VCF_merge/vqsr/recalibrated_snps_'+str(tranch)+'.vcf.gz'
	command = gatk_exe+' -T ApplyRecalibration -R '+ref_fasta+' -input '+vcf+' -mode SNP --ts_filter_level '+str(tranch)+' -recalFile '+recal_file+' -tranchesFile '+tranch_file+' -o '+output
	#print(gatk_exe+' -T ApplyRecalibration -R '+ref_fasta+' -input '+vcf+' -mode SNP --ts_filter_level '+str(tranch)+' -recalFile '+recal_file+' -tranchesFile '+tranch_file+' -o '+output)
	tranch = str(tranch)
	os.system('qsub -q evolgen.q -P evolgen -l rmem=20G -l mem=20G -b y -N rec'+tranch+' -e /fastdata/bo1tg/err/zf_ar'+str(tranch)+'.e -o /fastdata/bo1tg/out/zf_ar'+str(tranch)+'.o -cwd '+command)


