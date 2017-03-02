#!/usr/bin/env python
#code that generates the matrix.csv from the .var files

import sys
import collections
from collections import defaultdict
import os
import re

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

variants=set()
for line in open(sys.argv[1]):
    variants=line.rstrip().split(",")

d = defaultdict(set)
strains = set()
filename = sys.argv[2]
if filename.endswith(".var"):
    sn = filename.rsplit(".", 1)[0] #strain name
    for line in open(filename):
        e=line.rstrip().split("\t")
        d[sn].add(e[3]) #keys are strainnames and values are variantnames ie. snpnames
        strains.add(sn)
        #print >>sys.stderr, d[e[3]]

variants=list(set(variants))
strains=list(strains)
#print >>sys.stderr, sys.argv[1]

output2=sys.stdout
output2.write('strain,')
output2.write(",".join(variants))
output2.write("\n")
for s in strains:
 output2.write(s)
 for v in variants:
   parts=v.split('_')
   if parts[1] in ['CN','CS','CZ']: #coding snps pool changes that cause the same aachange
       #print parts[4]+"_"+parts[5]
       if re.search('\.$', parts[4]): #for to stop changes
	     parts[4]=parts[4][0:(len(parts[4])-1)] +"*"
       if re.search('^\.', parts[4]): #for stop to changes
	     parts[4]="*"+parts[4][1:]
       if re.search('\.$', parts[5]): #for oxyR'
             parts[5]=parts[5][0:(len(parts[5])-1)] +"'"	
       pattern=re.escape(parts[4])+"_"+re.escape(parts[5])	
       if re.search(pattern,"\t".join(d[s])):
             output2.write(",1")
             #print >>sys.stderr, v
       else:
             output2.write(",0")
   elif parts[1] in ['CF']: #coding frameshifts pool all that occur at the same nucleotide start
       pattern=re.compile(parts[2] + '_[^\s\,]+_' + parts[5])
       if re.search(pattern, "\t".join(d[s])):
             output2.write(",1")
#             print 'here!', s, v #>>sys.stderr, v
       else:
             output2.write(",0")
#	     print >>sys.stderr, v
   elif parts[1] == 'P': #promoter (to maintain compatibility with old naming used in randomforest built from MIP data
       operon=parts[5].split('.')
       pattern=parts[3]+'_'+parts[4]+'_'+operon[0]	
       if re.search(pattern, "\t".join(d[s])):
             output2.write(",1")
       else:
             output2.write(",0")
   else: #non coding, ribosomal changes SNPs or indels, and coding non frame shift indels have to be identical
       if v in d[s]:
           output2.write(",1")
           #print >>sys.stderr, v
       else:
           output2.write(",0")
 output2.write("\n")


