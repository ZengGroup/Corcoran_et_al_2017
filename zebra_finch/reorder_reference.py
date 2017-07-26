
import sys
sys.path.append('/home/bo1pgc/local/python')
import egglib

chrom_order = ['chr1', 'chr1A', 'chr1B', 'chr2', 'chr3', 'chr4', 'chr4A', 
               'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 
               'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 
               'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 
               'chr24', 'chr25', 'chr26', 'chr27', 'chr28', 'chrLG2', 
               'chrLG5', 'chrLGE22', 'chrZ', 'chrM', 'chrUn', 'chr1A_random', 
               'chr1B_random', 'chr1_random', 'chr2_random', 'chr3_random', 
               'chr4A_random', 'chr4_random', 'chr5_random', 'chr6_random', 
               'chr7_random', 'chr8_random', 'chr9_random', 'chr10_random', 
               'chr11_random', 'chr12_random', 'chr13_random', 'chr14_random', 
               'chr15_random', 'chr16_random', 'chr17_random', 'chr18_random', 
               'chr19_random', 'chr20_random', 'chr21_random', 'chr22_random', 
               'chr23_random', 'chr24_random', 'chr25_random', 'chr26_random', 
               'chr27_random', 'chr28_random', 'chrLGE22_random', 'chrZ_random']


ref = egglib.Container(sys.argv[1])

ref_order = ref.names()

reordered_ref = egglib.Container()

for i in chrom_order:
    index = ref_order.index(i)
    reordered_ref.append(ref[index][0], ref[index][1])
    
print reordered_ref




