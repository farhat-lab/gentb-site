#!/usr/bin/env python
"""
generates a matrix.csv file from one .var file
"""

import re
import os
import sys
import collections
from collections import defaultdict

from argparse import ArgumentParser, ArgumentTypeError

def file_type(fn):
    if os.path.isfile(fn):
        return fn
    raise ArgumentTypeError("File not found: %s" % fn)

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('mutations', help='CSV File of mutation names', type=file_type)
    parser.add_argument('strain', help='VAR File for a strain to process', type=file_type)
    args = parser.parse_args()

    with open(args.mutations, 'r') as fhl:
        variants = list(set(fhl.read().rstrip().split(",")))

    filename = args.strain

    d = {} #defaultdict(set)
    #strains = set()
    output2 = sys.stdout
    output2.write('strain,')
    output2.write(",".join(variants))
    output2.write("\n")

    sn = filename.rsplit(".", 1)[0] #strain name
    output2.write(sn)
    for line in open(filename):  
      a=line.rstrip().split("\t") #variant name i.e. snpname
      e=a[3]
      if e != "varname":
        parts=e.split('_')
        if e in variants:
            d[e]="1"
        elif parts[1] in ['CN','CS','CZ']: #coding snps pool changes that cause the same aachange
            if re.search('\*$', parts[4]): #for to stop changes
                parts[4]=parts[4][0:(len(parts[4])-1)] +"."
            if re.search('^\*', parts[4]): #for stop to changes
                parts[4]="."+parts[4][1:]
            if re.search("\'$", parts[5]): #for oxyR'
                parts[5]=parts[5][0:(len(parts[5])-1)] +"."
            pattern=parts[4]+"_"+parts[5]
            it = tuple(re.finditer(r"\s?([\w\.]+%s[\w\.]*)\s?"%pattern, "\t".join(variants), re.IGNORECASE))
            if len(it) >1 :
                d[it[0].group()]="1"
                #plan to add more matching to the right mutation here in a future version
            elif len(it)==1:             
                d[it[0].group()]="1"
        elif parts[1] in ['CF']: #coding frameshifts pool all that occur at the same nucleotide start
            pattern=parts[1] + '_' + parts[2] + '_[^\s\,]+_' + parts[5]
            it = tuple(re.finditer(r"\s?([\w\.]+%s[\w\.]*)\s?"%pattern, "\t".join(variants), re.IGNORECASE))
            if len(it) >1 :
                d[it[0].group()]="1"
                #plan to add more matching to the right mutation here in a future version
            elif len(it)==1:
                d[it[0].group()]="1"
        elif parts[1] == 'P': # promoter (to maintain compatibility with old naming used in randomforest built from MIP data
            #output2.write(parts[1] + "\n")
            operon=parts[len(parts)-1].split('-')
            if operon[0]=="promoter":
                pattern=parts[3]+'_'+operon[0]+'_'+operon[1]
            else:
                pattern=parts[3]+'_'+parts[4]+'_'+operon[0]
                it = tuple(re.finditer(r"\s?([\w\.]+%s[\w\.]*)\s?"%pattern, "\t".join(variants), re.IGNORECASE))
                if len(it) >1 :
                    d[it[0].group()]="1"
                    #plan to add more matching to the right mutation here in a future version
                elif len(it)==1:
                    d[it[0].group()]="1"

    for variant in variants:
        output2.write("," + str(int(variant in d)))
    output2.write("\n")

if __name__ == '__main__':
    main()

