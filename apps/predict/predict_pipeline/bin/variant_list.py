#code that generates the variant_name_list.csv

import os
import numpy as np
import random
import matplotlib.pyplot as plt
#%matplotlib inline
import collections
from collections import defaultdict

current=os.getcwd()
variants=set()
for filename in os.listdir(current):
    if filename.endswith(".var"):
        for line in open(filename):
                variants.add(line.rstrip().split("\t")[4])
                
output= file( "variant_name_list.csv", "w" )
output.write( ",".join(variants) )
output.write( "\n" )

#code that generates the matrix.csv for the mysql variants

variants=set()
for line in open("variant_name_list.csv"):
    variants=line.rstrip().split(",")

d = defaultdict(set)
strains=set()
for filename in os.listdir(current):
    if filename.endswith(".var"):
        for line in open(filename):
            e=line.rstrip().split("\t")
            d[e[4]].add(e[0]) #keys are variants names with values strain names
            strains.add(e[0])
#print(strains, file=sys.stderr)

variants=list(variants)
strains=list(strains)
output2=file('matrix.csv','w')
output2.write('strain,')
output2.write(",".join(variants))
output2.write("\n")
for s in strains:
 output2.write(s)
 for v in variants:
   if s in d[v]:
     output2.write(",1")
   else:
     output2.write(",0")
 output2.write("\n")
    

#code that generates the matrix.csv from the .var files

variants=set()
for line in open("variant_name_list.csv"):
    variants=line.rstrip().split(",")

d = defaultdict(set)
strains=set()
for filename in os.listdir(current):
    if filename.endswith(".var"):
        sn=filename.split(".")[0] #strain name
        for line in open(filename):
            e=line.rstrip().split("\t")
            d[e[4]].add(sn) #keys are variants names with values strain names
            strains.add(sn)
#print(strains, file=sys.stderr)

variants=list(variants)
strains=list(strains)
output2=file('matrix.csv','w')
output2.write('strain,')
output2.write(",".join(variants))
output2.write("\n")
for s in strains:
 output2.write(s)
 for v in variants:
   if s in d[v]:
     output2.write(",1")
   else:
     output2.write(",0")
 output2.write("\n")
    


            

