#!/usr/bin/env python

import sys
import os
from optparse import OptionParser

#python concatenate_gvcfs.py -v vcf -c $2 -l $QSUB_PATH -R $REF_FILE -p $VCF_FILE -n $3

parser = OptionParser()
parser.add_option('-v', help='Specify the type of VCF file, vcf or gvcf', dest='variant_file')
parser.add_option('-c', help='File listing chromosomes to be concatenated', dest='source_vcf')
parser.add_option('-p', help='Path where the concatenated VCF file will be written', \
                      dest='vcf_path')
parser.add_option('-n', help="name to append to the output concatenated vcf file", dest='sample_name')
parser.add_option('-l', help='Path to where the qsub log files will be written', \
                      dest='log_path')
parser.add_option('-R', help="Path to the reference genome fasta file", dest="ref_file")    
parser.add_option('-q', help="iceberg queue to submit to: general or evolgen", default="general", dest="iceberg_queue")
parser.add_option('-o', help="Output file name", dest="output_name")

(opts, args) = parser.parse_args()

if opts.log_path != '/':
    log_path = opts.log_path + '/'
else:
    log_path = opts.log_path

if opts.vcf_path[-1] != '/':
    vcf_path = opts.vcf_path + '/'
else:
    vcf_path = opts.vcf_path
    


chrom_order = [] 
for i in open(opts.source_vcf): # use a file listing chromosomes
    chrom = i.strip()
    if chrom == 'scaffolds.intervals':
        chrom_order.append('scaffolds')    
    else:
        chrom_order.append(chrom)


reference_file= "-R " + opts.ref_file + ' '

variant_files = ''

sample = opts.sample_name
sample_outfile = opts.output_name


for chrom in chrom_order:
    if opts.variant_file == 'gvcf':
        variant_string = (' -V ' + vcf_path + sample + '.' + chrom  +  '.raw.snps.indels.g.vcf')

    elif opts.variant_file == 'vcf':
        variant_string = (' -V ' + vcf_path + sample + '.' + chrom  +  '.raw.snps.indels.vcf')
        
    else:
        print 'Error: Need to specify vcf or gvcf'
        sys.exit(1)

    variant_files += variant_string



if opts.variant_file == 'gvcf':
    variant_string = (' -V ' + vcf_path + sample + '.' + chrom  +  '.raw.snps.indels.g.vcf')
    combined_vcf = ' -out ' + vcf_path + sample_outfile
    rm_pattern_1 = vcf_path + sample + '.*.raw.snps.indels.g.vcf'
    rm_pattern_2 = vcf_path + sample + '.*.raw.snps.indels.g.vcf.idx'
    combine_vcf = ("java -Xmx8g -cp ~/bin/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants "
                   + reference_file + variant_files + combined_vcf + " -assumeSorted")
    hold_id = "-hold_jid hap_call." + sample + '.* '


elif opts.variant_file == 'vcf':
    variant_string = (' -V ' + vcf_path + sample + '.' + chrom  +  '.raw.snps.indels.vcf')
    combined_vcf = ' -out ' + vcf_path + sample_outfile
    rm_pattern_1 = vcf_path + sample + '.*.raw.snps.indels.vcf'
    rm_pattern_2 = vcf_path + sample + '.*.raw.snps.indels.vcf.idx'
    combine_vcf = ("java -Xmx8g -cp ~/bin/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants "
                   + reference_file + variant_files + combined_vcf + " -assumeSorted")
    hold_id = "-hold_jid genotype." + sample + '.* '    



rm_part_vcf = "rm " + rm_pattern_1 + ' ' + rm_pattern_2
    
    
gather_job_id = "-N gather_vcf." + sample
    
e_log = "-e " + log_path + 'gather_vcf.' + sample + '.e '
o_log = "-o " + log_path + 'gather_vcf.' + sample + '.o '

if opts.iceberg_queue == 'general':
    queue = ""
elif opts.iceberg_queue == 'evolgen':
    queue = "-P evolgen -q evolgen.q "
    
qsub_gather = ("qsub -b y " + queue + gather_job_id + 
               " -l mem=12G -l rmem=12G  -l arch=intel* -l h_rt=24:00:00 " + 
               hold_id + e_log + o_log + "-cwd " + combine_vcf)

os.system(qsub_gather)

rm_job_id = "-N rm_part_vcf." + sample
rm_hold_id = "-hold_jid gather_vcf." + sample
rm_e_log = " -e " + log_path + "rm_part_vcf." + sample + '.e '
rm_o_log = "-o " + log_path + "rm_part_vcf." + sample + '.o '


qsub_rm = ("qsub -b y " + queue + rm_job_id +
           " -l mem=12G -l rmem=12G  -l arch=intel* -l h_rt=8:00:00 " +
           rm_hold_id + rm_e_log + rm_o_log + "-cwd " + rm_part_vcf)

#os.system(qsub_rm) 
#print qsub_gather
#print qsub_rm
