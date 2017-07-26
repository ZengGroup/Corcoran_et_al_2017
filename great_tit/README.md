# Steps in the processing of the greet tit BGI sequence data

This document describes the mapping and SNP calling from paired end Illumina sequence data from 10 great 
tit individuals sample across Europe. Each bird was seqenced on a single lane of HiSeq 2500 machine, resulting in  ~60Gb of sequence data per individual. 


## Public Availibility of sequence data and files
The fastq file used are publically available at [here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA381923)

The Version of the great tit genome (V1.04) used in this pipeline is publically available [here](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/522/545/GCF_001522545.1_Parus_major1.0.3)

## Programs used in this pipeline
__Versions of programs used__

* cutadapt v1.8.1
* FastQC 0.11.3 
* samtools 1.2
* bwa 0.7.12-r1039
* Picard 1.130
* GATK 3.4-46-gbc02625
* bcftools 1.2
* freebayes 1.02
* faCount and faList from https://www.biostat.wisc.edu/~cdewey/software/cndsrc-2013.01.11.tar.gz
* vcftools v0.1.14

## Steps in the pipeline

### 1. Sequence data preparation and reference genome files

#### Preparing the reference genome

The reference genome used was indexed in the same way as for the usng bwa
as follows:

```
$ sh prepare_reference.sh
```

The cleaned fastq files received from the BGI were copied to the following 
location /fastdata/bo1pgc/parus_reseq/bgi_reseq/fastq_files. The names of the
fastq files for each sample are listed in Table 2 below.

_Table 1_ Cleaned fastq files provided by BGI and sample ids for the birds

| Fastq_1 | Fastq_2 | Sample | Sex |
|:--|:--|:--|:--:|
|FCC7GKLANXX_L8_wHAIPI021608-93_1.fq.gz | FCC7GKLANXX_L8_wHAIPI021608-93_2.fq.gz | 1280 | Male |
|FCC7GT8ANXX_L6_wHAIPI021614-95_1.fq.gz | FCC7GT8ANXX_L6_wHAIPI021614-95_2.fq.gz | 1485 | Male |
|FCC7GKLANXX_L6_wHAIPI021616-90_1.fq.gz | FCC7GKLANXX_L6_wHAIPI021616-90_2.fq.gz | 15 | Male |
|FCC7FYRANXX_L5_wHAIPI021634-89_1.fq.gz | FCC7FYRANXX_L5_wHAIPI021634-89_2.fq.gz | 167 | Male |
|FCC7GT8ANXX_L8_wHAIPI021636-97_1.fq.gz | FCC7GT8ANXX_L8_wHAIPI021636-97_2.fq.gz | 249 | Male |
|FCC7FYRANXX_L4_wHAIPI021611-88_1.fq.gz | FCC7FYRANXX_L4_wHAIPI021611-88_2.fq.gz | 318 | Male |
|FCC7GKLANXX_L7_wHAIPI021620-92_1.fq.gz | FCC7GKLANXX_L7_wHAIPI021620-92_2.fq.gz | 61 | Male |
|FCC5LDHANXX-PARxfiRAADIAAPEI-42_L3_1.fq.gz | FCC5LDHANXX-PARxfiRAADIAAPEI-42_L3_2.fq.gz | 917 | Male |
|FCC7GT8ANXX_L7_wHAIPI021624-96_1.fq.gz | FCC7GT8ANXX_L7_wHAIPI021624-96_2.fq.gz | 943 | Male |
|FCC7GT8ANXX_L5_wHAIPI021613-94_1.fq.gz | FCC7GT8ANXX_L5_wHAIPI021613-94_2.fq.gz | TR43666 | Male |


Each sample has a folder with them name of the sample (e.g. 917) and within each
folder there is a folder called clean_reads which contains the clean reads in 
compressed fastq format. The raw reads are found in the folder called raw_reads. 
Both the raw and cleaned read are permanently stored on the evolgen1 share_data 
storage, accessible from iceberg.

#### Constructing the read group information
The read group information for each pair of fastq file from the same lane, 
library and sample was constructed using the ```create_fq_and_rg_bgi.sh``` 
shell script. This script requires the path to the fastq files to be given within 
the script and will produce a tab-delimited file with the fastq files and associated 
read group information, needed by downstream tools, on each row. The paths to fastq
files needs to be specified in the shell script. This script produces two text files.
The first called fastq_and_rg.txt and the second called fastq_and_rg_clean.txt


    $ sh create_fq_and_rg_bgi.sh

#### Checking the quality of the reads

The quality of the reads in the fastq files were checked with the FastQC program


    $ sh fastqc.sh fastq_and_rg.txt


The plots of kmers indicated that some kmers were over-represented at 
the beginning of the reads and represented some remaining adapter contamination. 
These adapter sequences were removed using the cutadapt program called 
in the script called ``` filter_adapters.sh```


    $ sh filter_adapters.sh fastq_and_rg.txt


### 2. Mapping reads

The cleaned reads in the gzipped fastq files produced by BGI were mapped to 
the great tit reference genome (Version 1.04) for each lane of sequencing 
ssing the ```bwa_mem_real.sh``` script. The ```bwa_mem_real.sh``` script maps 
the read with bwa mem. It then performs the SAM to BAM conversions, BAM sorting 
and the removal of duplicates using Picard tools. The coverage of the sequencing
is checked using WgsMetrics from Picard. Plots and tables are written to a folder
called plots using the ```get_depth_stats.sh```. The R script used to generate the
plots requires ggplot2 to be installed for the relevant version of R that is loaded
as a module on iceberg. A short markdown report containing a table and the plots was 
also generated.


    $ sh bwa_mem_real.sh fastq_and_rg_clean.txt


The realigned BAM file were run through the picard tool CollectInsertSizeMetrics 
to get summary statistics on the insert sizes and to generate a histogram 
of the insert size distribution


    $ ls /fastdata/bo1pgc/parus_reseq/bgi_reseq/bams/*.bam > unrecalibrated_bams_bgi.txt
    $ sh collect_insert_size.sh unrecalibrated_bams_bgi.txt

#### Summary of Coverage

The summary of the Coverage across the genome was gathered using the ```get_depth_stats.sh``` script. The path to the BAM files needs to be specified inside the script.

    $ qsub get_depth_stats.sh

#### Sex determination

The 10 birds sequenced are supposed to be male, but I checked for any females, as they need to be treated differently when calling varinats and genotyping the z chromosome, as the haplotype caller need to be run in haploid mode for the Z chromosome in females. To identity any females, the realigned BAM files produced in previous step were used to calculatethe mean coverage on autosomes and Z chromosome for each bird using the samtools depth command. I then calculated the ratio of mean Z depth to Autosome depth. An ratio of ~0.5 would be expected ratio for a female sample. To do this I first generated a file called genome_file.txt that listed the chromosome name, chromsome length and mumber of non-N bases on the chromosome on each row of the file using the faCount program available from [here](https://www.biostat.wisc.edu/~cdewey/software.html)


    $ faCount < /data/bo1pgc/gtit_1.04/genome/Parus_major_1.04.rename.fa | while read i; do grep -v ^# | awk '{print $1, $2, $2-$7}'; done | grep -v total > genome_file.txt
    $ sh calc_Z_to_A_coverage.sh unrecalibrated_bams_wash.txt genome_file.txt z_to_aut_cov.txt


All samples showed no indication of being females, in concordance with the sexing of the same samples using the SNP array data on the Z chromosome for the same 10 samples. 
This information was added to the ```unrecalibrated_bams_bgi.txt```  and the resulting file has the following tab delimited format:

```
/fastdata/bo1pgc/parus_reseq/bgi_reseq/bams/1280.dedup.real.bam	M
/fastdata/bo1pgc/parus_reseq/bgi_reseq/bams/1485.dedup.real.bam	M
/fastdata/bo1pgc/parus_reseq/bgi_reseq/bams/15.dedup.real.bam	M
/fastdata/bo1pgc/parus_reseq/bgi_reseq/bams/167.dedup.real.bam	M
/fastdata/bo1pgc/parus_reseq/bgi_reseq/bams/249-R.dedup.real.bam	M
/fastdata/bo1pgc/parus_reseq/bgi_reseq/bams/318.dedup.real.bam		M
/fastdata/bo1pgc/parus_reseq/bgi_reseq/bams/61.dedup.real.bam		M
/fastdata/bo1pgc/parus_reseq/bgi_reseq/bams/917.dedup.real.bam		M
/fastdata/bo1pgc/parus_reseq/bgi_reseq/bams/943-R.dedup.real.bam	M
/fastdata/bo1pgc/parus_reseq/bgi_reseq/bams/TR43666.dedup.real.bam	M
```

### 3. Base Quality Score Recalibration (BQSR)

The realigned BAM file were run through the GATK HaplotypeCaller to generate gVCFs for each individual. The gVCF files were run through the GenotypeGvcf GATK tool to generate a multisample raw VCF file containing SNPs only. The HaplotypeCaller is run on each chromosome and the chromsomal gVCF files were conctaenated with the concatenate_vcf.py script which calls the GATK CatVariants tool. Chromosomes are listed in chromosome.txt and scaffolds in scaffold.interval


    $ faList < /data/bo1pgc/gtit_1.04/genome/Parus_major_1.04.rename.fa | grep ^Scaffold > scaffold.interval
    $ sh haplotype_caller.sh unrecalibrated_bams_bgi.txt chromosomes.txt
    $ ls /fastdata/bo1pgc/parus_reseq/bgi_reseq/unrecal_gvcf/*.vcf.gz > unrecalibrated_gvcf.bgi.list
    $ sh genotype_gvcfs.sh unrecalibrated_gvcf.bgi.list chromosomes.txt bgi_10birds


This raw VCF fle was then hard filtered using the VarianFiltration tool in GATK within the hard_filter_vcf.sh script. The hard filters used are specified within the script and follow recommendations on the GATK online documentation found [here](http://gatkforums.broadinstitute.org/discussion/2806/howto-apply-hard-filters-to-a-call-set)


    $ sh hard_filter_vcf.sh /fastdata/bo1pgc/parus_reseq/bgi_reseq/unrecal_vcf/bgi_10birds.raw.snps.indels.vcf.gz


The BQSR step of the GATK was carried out using the hard filtered SNP and INDEL vcf files, generated from the hard_filter_vcf.sh script above, as the know variant site. The plots of base quality scores before and after recalibration were generated using the AnalyzeCovariates tool in GATK.


    $ sh bqsr.sh unrecalibrated_bams_bgi.txt /fastdata/bo1pgc/parus_reseq/bgi_reseq/vcfs_for_bqsr/bgi_10birds.gatk_hard_filter.snp.vcf.gz /fastdata/bo1pgc/parus_reseq/bgi_reseq/vcfs_for_bqsr/bgi_10birds.gatk_hard_filter.indel.vcf.gz
    $ qsub analyse_covariates.sh unrecalibrated_bams_bgi.txt


#### Summary of Coverage

The summary of the Coverage after BQSR was gathered using the ```get_depth_stats.sh``` script. The path to the BAM files needs to be specified inside the script.

    $ sh calc_wgs_metrics.sh
    $ qsub get_depth_stats.sh
    $ sh get_mapping_stats.sh recalibrated_bams_bgi.txt


### 4. Variant Calling and genotyping  with GATK

The recalibrated BAM files for the 10 Birds were used to generate vcf files contain variants and genotype calls using the HaplotypeCaller and GenotypeGVCFs tools in GATK. The sex information on each bird was also included in the ```recalibrated_bams_bgi.txt``` file generated with `ls` command below:


    $ ls /fastdata/bo1pgc/parus_reseq/bgi_reseq/recalibrated_bams/*.bam > recalibrated_bams_bgi.txt
    $ sh haplotype_caller.sh recalibrated_bams_bgi.txt chromosomes.txt


The genotyping was perforemd using the GenotypeGVCFs tool. The allSites VCF file produced had a number of monomorphic (i.e., non-variant) sites that had an ALT allele present. These unwanted ALT alleles were removed using the SelectVariants tool, as these sites were being considered as SNPs in VQSR eventhough they are not SNPs. These sites were not emitted when a variant only VCF was produced. Apart from this the variant only and all-sites VCF produced from the same gVCF files had the same number of indels and SNPs output. 


    $ ls /fastdata/bo1pgc/parus_reseq/bgi_reseq/recal_gvcf/*.vcf.gz > recalibrated_gvcf.bgi.list
    $ sh genotype_gvcfs.sh recalibrated_gvcf.bgi.list chromosomes.txt bgi_10birds


### 5. Variant Quality Score Recalibration

#### Generating the training set

Called SNPs with Freebayes using the following command:
    
    $ freebayes+' -f  -L recalibrated_bams_bgi.txt  --report-monomorphic --report-genotype-likelihood-max --min-mapping-quality 20 --min-base-quality 10
    
Extracted and filtered SNPs based on depth and Quality

	$ python extract_SNPS_only.py
	$ qsub sort_and_bgzip.sh

Intersected the filtered GATK SNPs with the filtered Freebayes SNPs to generate the training set using bcftools isec.

    $ bcftools isec -p /fastdata/bo1tg/VCF_10birds/VCF_merge/Merged_VCF /fastdata/bo1tg/VCF_10birds/VCF_merge/gt_10birds.gatk.SNPsONLY.sorted.vcf.gz /fastdata/bo1tg/VCF_10birds/VCF_10birds/VCF_merge/gt_10birds.freebayes.SNPsONLY.sorted.vcf.gz -n=2 -w 1 
    
#### Running and applying VQSR to the GATK callset
Used training set to run VQSR and applied the recalibration tables generated by the VariantRecalibrator tool in GATK on the
unfiltered GATK callset. 

    $ qsub vqsr.sh
    $ python apply_recalibration.py
    
    
### 6. Repetitive Region filter

Marked sites in the recalibrated VCF files that were located in repetitive regions of the genome.

    $ qsub filter_repeats.gt.sh /fastdata/bo1tg/VCF_Backup_Padraic/great_tit/All_sites/recalibrated_snps_99.9.vcf.gz /fastdata/bo1tg/VCF_Backup_Padraic/great_tit/All_sites/recalibrated_snps_99.9.rep_filtered.vcf.gz 
        
