#code that generates the matrix.csv from the .var files

import sys
import collections
from collections import defaultdict
import os

variants=set()
for line in open("/groups/murray/run_pipeline/bin/variant_name_list.csv"):
    variants=line.rstrip().split(",")

d = defaultdict(set)
strains=set()
for filename in os.listdir(sys.argv[1]):
#    print >>sys.stderr, filename
    if filename.endswith(".var"):
        sn=filename.split(".")[0] #strain name
	filen=sys.argv[1]+"/"+filename
        for line in open(filen):
            e=line.rstrip().split("\t")
            d[sn].add(e[3]) #keys are strainnames and values are variantnames ie. snpnames
            strains.add(sn)
	    #print >>sys.stderr, d[e[3]]

variants=list(variants)
strains=list(strains)
#print >>sys.stderr, sys.argv[1]

f=sys.argv[1]+"/matrix.csv"
output2=file(f,'w')
output2.write('strain,')
output2.write(",".join(variants))
output2.write("\n")
for s in strains:
 #print >>sys.stderr, s
 output2.write(s)
 for v in variants:
   #print >>sys.stderr, v
   if v in d[s]:
     #print >>sys.stderr,d[s]
     output2.write(",1")
   else:
     output2.write(",0")
 output2.write("\n")
    


            

