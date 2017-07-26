import sys
import time
import gzip

fq_path = sys.argv[1]
fq_clean_path = sys.argv[2]

if fq_path[-1] != '/':
       fq_path += '/'

if fq_clean_path[-1] != '/':
       fq_clean_path += '/'


#generate time for dt read group tag
#should be date of sequencing, but just use current
# date to fill in this tag, so that the file produced
# works with the bwa_mem_real.sh script for mapping
date = time.localtime() 
dt = str(date[0]) + '-' + str(date[1]) + '-' + str(date[2]) 
    
f1 = open('fastq_and_rg.txt', 'w')
f2 = open('fastq_and_rg_clean.txt', 'w')
for i in sys.stdin:
    
        
    col = i.strip().split('/')
    
    sm = col[0]
    if  col[1][-7] == '1': 
        fq_id = col[1][0:-8]
    else:
        continue


    #get the flowcell and lane id from the first record of fastq 1
    fq_1 = fq_path + sm + '/clean_reads/' + fq_id + '_1.fq.gz'     
    fq_2 = fq_path + sm + '/clean_reads/' + fq_id + '_2.fq.gz'
    
    fq = gzip.open(fq_1, 'rb')
    fq_header = fq.readline()
    fq.close()
    
    fq_fields = fq_header.split(':')
        
    flowcell = fq_fields[2]
    lane = fq_fields[3]
    
    pl = 'illumina'
    pu = flowcell + '_' + lane
    lb = sm + '_' + 'lib1'
    pi = '500'
    ds = 'paired_end'
    cn = 'BGI'

    rg_id = lb + '_' + lane

    #print fq_1 + '\t' + fq_2 + '\t' + rg_id + '\t' + pl \
    #+ '\t' + pu + '\t' + lb + '\t' + pi + '\t' + \
    #ds + '\t' + dt + '\t' + sm + '\t' + cn
        
        
    f1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % 
            (fq_1, fq_2, rg_id, pl, pu, lb, pi, ds, dt, sm, cn))

    
    fq_clean_1 = fq_clean_path + sm  + '.clean' + '.1.fq.gz'
    fq_clean_2 = fq_clean_path + sm  + '.clean' + '.2.fq.gz'
    
    f2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
          (fq_clean_1, fq_clean_2, rg_id, pl, pu, lb, pi, ds, dt, sm, cn))
