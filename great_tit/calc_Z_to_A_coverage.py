#!/usr/bin/env python

import sys


sample = sys.argv[1].split('/')[-1].split('.')[0]

aut_depth = 0
z_depth = 0

aut_bases = 0
z_bases = 0

aut_done = 0
# get the number of non-N sites in the genome for autosomes and Z

for i in open(sys.argv[2]):
    col = i.rstrip().split()
    
    if col[0].startswith('Scaffold'):
        continue
    elif col[0] != 'chrZ':
        aut_bases += int(col[2])
    elif col[0] == 'chrZ':
        z_bases += int(col[2])
    
for line in sys.stdin:
    col = line.rstrip().split()    
    if col[0].startswith('Scaffold'):
        continue
    elif col[0] != 'chrZ':
        aut_depth += int(col[2])  
    elif col[0] == 'chrZ':
        z_depth += int(col[2])


mean_aut = aut_depth / float(aut_bases)
mean_z = z_depth / float(z_bases)

z_aut_ratio =  mean_z / mean_aut

print sample + '\t' + str(mean_aut) + '\t' \
+ str(mean_z) + '\t' + str(z_aut_ratio)
