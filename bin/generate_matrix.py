#!/usr/bin/env python
#code that generates a matrix.csv file from one .var file

import sys
import collections
from collections import defaultdict
import os
import re

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

variants=set()
for line in open(sys.argv[1]):
    variants=line.rstrip().split(",")

variants=list(set(variants))

d = {} #defaultdict(set)
#strains = set()
output2=sys.stdout
output2.write('strain,')
output2.write(",".join(variants))
output2.write("\n")

filename = sys.argv[2]
if filename.endswith(".var"):
    sn = filename.rsplit(".", 1)[0] #strain name
    output2.write(sn)
    for line in open(filename):  
      a=line.rstrip().split("\t") #variant name i.e. snpname
      e=a[3]
      if e != "varname":
        #sys.stderr.write(str(e)+"\n")
        parts=e.split('_')
        if e in variants:
            d[e]="1"
	    #sys.stderr.write("match\n")
        elif parts[1] in ['CN','CS','CZ']: #coding snps pool changes that cause the same aachange
            if re.search('\*$', parts[4]): #for to stop changes
                parts[4]=parts[4][0:(len(parts[4])-1)] +"."
            if re.search('^\*', parts[4]): #for stop to changes
                parts[4]="."+parts[4][1:]
            if re.search("\'$", parts[5]): #for oxyR'
                parts[5]=parts[5][0:(len(parts[5])-1)] +"."
            pattern=parts[4]+"_"+parts[5]
            it = tuple(re.finditer(r"\s?([\w\.]+%s[\w\.]*)\s?"%pattern, "\t".join(variants), re.IGNORECASE))
	    #sys.stderr.write(it[0].group())
	    if len(it) >1 :
		d[it[0].group()]="1"
		#plan to add more matching to the right mutation here in a future version
#		sys.stderr.write("match1+\n")
	    elif len(it)==1:     	
                d[it[0].group()]="1"
#		sys.stderr.write("match1\n")
        elif parts[1] in ['CF']: #coding frameshifts pool all that occur at the same nucleotide start
            pattern=parts[1] + '_' + parts[2] + '_[^\s\,]+_' + parts[5]
	    it = tuple(re.finditer(r"\s?([\w\.]+%s[\w\.]*)\s?"%pattern, "\t".join(variants), re.IGNORECASE))
            if len(it) >1 :
                d[it[0].group()]="1"
                #plan to add more matching to the right mutation here in a future version
#		sys.stderr.write("match2+\n")
            elif len(it)==1:
                d[it[0].group()]="1"
#	        sys.stderr.write('match2'+"\n")
        elif parts[1] == 'P' and False: #promoter (to maintain compatibility with old naming used in randomforest built from MIP data
	   # sys.stderr.write(parts[2] + "\n")
            operon=parts[5].split('-')
            pattern=parts[3]+'_'+parts[4]+'_'+operon[0]
            it = tuple(re.finditer(r"\s?([\w\.]+%s[\w\.]*)\s?"%pattern, "\t".join(variants), re.IGNORECASE))
            if len(it) >1 :
                d[it[0].group()]="1"
                #plan to add more matching to the right mutation here in a future version
#		sys.stderr.write("match3+\n")
            elif len(it)==1:
                d[it[0].group()]="1"
#	        sys.stderr.write("match3\n")
        
#sys.stderr.write("\t".join(d.keys()))

for v in variants:
    if v in d.keys():
        output2.write(",1")
    else:
        output2.write(",0")
output2.write("\n")


