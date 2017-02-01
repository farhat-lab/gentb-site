#code that generates the matrix.csv from the .var files

import sys
import collections
from collections import defaultdict
import os
import re

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
#SCRIPT_DIR = '/home/maha/gentb/gentb-site/apps/predict/pipeline/bin'

variants=set()
for line in open(SCRIPT_DIR + "/variant_name_list.csv"):
    variants=line.rstrip().split(",")

d = defaultdict(set)
strains=set()
for filename in os.listdir(sys.argv[1]):
#for filename in os.listdir('/home/maha/gentb/gentb-site/apps/predict/pipeline/bin'):
#    print >>sys.stderr, filename
    if filename.endswith(".var"):
        sn=filename.split(".")[0] #strain name
        filen=sys.argv[1]+"/"+filename
        #filen='/home/maha/gentb/gentb-site/apps/predict/pipeline/bin'+"/"+filename
        for line in open(filen):
            e=line.rstrip().split("\t")
            d[sn].add(e[3]) #keys are strainnames and values are variantnames ie. snpnames
            strains.add(sn)
            #print >>sys.stderr, d[e[3]]

variants=list(set(variants))
strains=list(strains)
#print >>sys.stderr, sys.argv[1]

f=sys.argv[1]+"/matrix.csv"
output2=file(f,'w')
output2.write('strain,')
output2.write(",".join(variants))
output2.write("\n")
for s in strains:
 output2.write(s)
 for v in variants:
   parts=v.split('_')
   if parts[1] in ['CN','CS','CZ']: #coding snps pool changes that cause the same aachange
       for u in d[s]:
           if re.search(parts[5],u) and re.search(parts[4],u):
             output2.write(",1")
           else:
             output2.write(",0")
   elif parts[1] in ['CF']: #coding frameshifts pool all that occur at the same nucleotide start
       for u in d[s]:
           if re.search(parts[5], u) and re.search(parts[2],u):
             output2.write(",1")
           else:
             output2.write(",0")
   else: #non coding, promoter, ribosomal changes SNPs or indels, and coding non frame shift indels have to be identical
       if v in d[s]:
           output2.write(",1")
       else:
           output2.write(",0")
 output2.write("\n")


