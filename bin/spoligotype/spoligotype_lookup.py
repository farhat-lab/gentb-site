#!/usr/bin/env python

# Copyright 2014, Thomas R. Ioerger ..
import sys,os,re

def octal(s): # converts to binary
  o = ''
  for x in s[:-1]:
    if x=='0': o += '000'
    if x=='1': o += '001'
    if x=='2': o += '010'
    if x=='3': o += '011'
    if x=='4': o += '100'
    if x=='5': o += '101'
    if x=='6': o += '110'
    if x=='7': o += '111'
  o += s[-1] # last bit, 43rd spacer
  return o

def subsumes(a,b):
  if len(a)!=43 or len(b)!=43: print "Error in subsumes"; sys.exit(0)
  good = True
  for i in range(len(a)):
    if a[i]=='0' and b[i]=='1': good = False
  return good

def similarity(a,b):
  count = 0
  for i in range(len(a)):
    if a[i]==b[i]: count += 1
  return count

if __name__=="__main__":

# detects the octal representation created upstream
  for line in sys.stdin:
    if line.startswith("! Octal:"):
      fname = re.findall('\[.*\]', line)[0]
      if not fname:
        sys.exit("Error: No file name detected.")
      query = line[8 + len(fname):].rstrip() # remove trailing newline char and fname
      fname = re.split('/', fname[1:-1])[-1][:-6]
      # remove braces, split string by /, take last entry (file name) and remove fastq portion
      break
  try:
    binquery = octal(query)
  except NameError:
    sys.exit("Error: No octal representation found in input.")
  
  lines = open(sys.argv[1]).readlines()
  lines = [x.rstrip() for x in lines]
  data = [x.split() for x in lines]
  
  scores = []
  for i in range(len(data)):
    id,spol = data[i][0],data[i][1]
    octspol = octal(spol)
    fam = ' '.join(data[i][2:])
    if subsumes(octspol,binquery)==True:
      sim = similarity(octspol,binquery)
      scores.append((sim,i,id,spol,octspol,fam))
  
  scores.sort()
  scores.reverse()
  
  mostsim = scores[:1][0]
  # mostsim is a tuple with order (sim,idx,id,spol,octspol,fam)
  output = [fname, query, mostsim[2], mostsim[3], mostsim[0], mostsim[5]]
  print('\t'.join(map(str,output)))
  
'''
  print "Query",query,binquery
  print
  print 'top SUBSUMING matches sorted by similarity:'
  for (sim,idx,id,spol,octspol,fam) in scores[:5]:
    print "%-5s %s %s sim=%d %s" % (id,spol,octspol,sim,fam)

  print
  print 'max similarity is 43 (identical)'
  print 'first number is SpolDb4 strain index'
  print 'end of line is strain family'
'''

''' Example output
File Name	Spoligotype	RefID	SpolDb4 strain index	Similarity (Max = 43)	Family
5	775776774020731	1742	775777774020771	41	H1
6	775776774020731	1742	775777774020771	41	H1
'''
