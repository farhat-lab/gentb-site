#!/usr/bin/env python
"""
generates a matrix.csv file from one .var file
"""

import re
import os
import sys

from argparse import ArgumentParser, ArgumentTypeError

def file_type(fname):
    """Confirm file input is in a file"""
    if os.path.isfile(fname):
        return fname
    raise ArgumentTypeError("File not found: %s" % fname)

def main():
    """Main body of the script"""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('mutations', help='CSV File of mutation names', type=file_type)
    parser.add_argument('strain', help='VAR File for a strain to process', type=file_type)
    args = parser.parse_args()

    with open(args.mutations, 'r') as fhl:
        variants = list(set(fhl.read().rstrip().split(",")))
    generate_matrix(variants, args.strain)

def generate_matrix(variants, filename):
    """Generate the matrix"""
    out = {}
    #strains = set()
    output2 = sys.stdout
    output2.write('strain,')
    output2.write(",".join(variants))
    output2.write("\n")
    output2.write(filename.rsplit(".", 1)[0]) # strain name

    for line in open(filename):
        items = line.rstrip().split("\t") #variant name i.e. snpname
        snpname = items[5]
        if snpname != "varname":
            parts = snpname.split('_')
            if snpname in variants:
                out[snpname] = "1"
                #sys.stderr.write(snpname)
                #sys.stderr.write("\tmatch!!!!!!!\n")
            elif parts[1] in ['CN', 'CS', 'CZ']:
                # coding snps pool changes that cause the same aachange
                if re.search(r'\*$', parts[4]): #for to stop changes
                    parts[4] = parts[4][0:(len(parts[4])-1)] + "."
                if re.search(r'^\*', parts[4]): #for stop to changes
                    parts[4] = "." + parts[4][1:]
                if re.search(r"'$", parts[5]): #for oxyR'
                    parts[5] = parts[5][0:(len(parts[5])-1)] + "."
                pattern = parts[4] + "_" + parts[5]
                pattern = r"\s?([\w\.]+%s[\w\.]*)\s?" % pattern.replace('*', r'\*')
                #print(pattern)
                ita = tuple(re.finditer(pattern, ",".join(variants))) #, re.IGNORECASE))
                if len(ita) > 1:
                    out[ita[0].group()] = "1"
                    #plan to add more matching to the right mutation here in a future version
                    #sys.stderr.write(ita[0].group())
                    #sys.stderr.write("match!!!!!!!\n")
                elif len(ita) == 1:
                    out[ita[0].group()] = "1"
                    #sys.stderr.write(ita[0].group())
                    #sys.stderr.write("match!!!!!!!\n")
            elif parts[1] in ['CF']:
                #coding frameshifts pool all that occur at the same nucleotide start
                pattern = parts[1] + '_' + parts[2] + r'_[^\s\,]+_' + parts[5]
                pattern = r"\s?([\w\.]+%s[\w\.]*)\s?" % pattern.replace('*', r'\*')
                itb = tuple(re.finditer(pattern, ",".join(variants))) #, re.IGNORECASE))
                if len(itb) > 1:
                    out[itb[0].group()] = "1"
                    #plan to add more matching to the right mutation here in a future version
                    #sys.stderr.write(itb[0].group())
                    #sys.stderr.write("match!!!!!!!\n")
                elif len(itb) == 1:
                    out[itb[0].group()] = "1"
                    #sys.stderr.write(itb[0].group())
                    #sys.stderr.write("match!!!!!!!\n")
                #if any other frameshift mutation is found that is not in the variant list
                #write a 1 in the respective field of another variant that is part of the variant list
                elif len(itb) == 0:
                    if parts[5] in ['pncA']:
                        out['INS_CF_2288725_i517C_173_pncA'] = "1"
                    
            elif parts[1] == 'P':
                # promoter (to maintain compatibility with old naming used in
                # randomforest built from MIP data
                operon = parts[-1].split('-')
                #sys.stderr.write(str(operon))
                if operon[0] == "promoter":
                    pattern = parts[3] + '_' + operon[0] + '_' + operon[1]
                else:
                    pattern = parts[3] + '_' + parts[4] + '_' + operon[0]

                #sys.stderr.write(str(pattern))
                #sys.stderr.write("\n")
                pattern = r"\s?([\w\.]+%s[\w\.]*)\s?" % pattern.replace('*', r'\*')
                itc = tuple(re.finditer(pattern, ",".join(variants))) #, re.IGNORECASE))
                if len(itc) > 1:
                    out[itc[0].group()] = "1"
                    #plan to add more matching to the right mutation here in a future version
                elif len(itc) == 1:
                    out[itc[0].group()] = "1"
                    #sys.stderr.write(itc[0].group())
                    #sys.stderr.write("match!!!!!!!\n")

            #include all variants 40 basepairs up or downstream of the rifampicin resistance determining region
            elif ( parts[1] in ['CF', 'CI', 'CD', 'CN', 'CZ'] ) and ( '761041' <= parts[2] <= '761202' ) and ( parts[5] in ['rpoB'] ):
                out['SNP_CN_761155_C1349T_S450L_rpoB'] = "1"

            #match non-synonymous SNPs in pncA to a variant that the model has seen 
            elif ( parts[1] in ['CN'] ) and ( parts[5] in ['pncA'] ):
                out['SNP_CN_2289090_T152C_H51R_pncA'] = "1"

            #match all deletions and stopcodon mutations to a stopcodon mutation that the model has seen
            elif ( parts[1] in ['CD, CZ'] ) and ( parts[5] in ['pncA'] ):
                out['SNP_CZ_2288933_G309C_Y103._pncA'] = "1" 

            #match all stopcodons variants in katG to a katG stopcodon that the model has seen
            elif ( parts[1] in ['CZ', 'CF'] )  and ( parts[5] in ['katG'] ):
                out['SNP_CZ_2155101_G1011C_Y337._katG'] = "1"

            #match all stopcodon and frameshift variants in ethA to a stopcodon variant the model has seen
            elif ( parts[1] in ['CF', 'CZ'] ) and ( parts[5] in ['ethA'] ):
                out['SNP_CZ_4326714_G760A_Q254._ethA'] = "1"

            #match all frameshifts in tlyA 
            elif ( parts[1] in ['CF'] ) and ( parts[5] in ['tlyA'] ):
                out['SNP_CN_1918160_T221A_L74Q_tlyA'] = "1"

            #match all frameshift mutations in gid
            elif ( parts[1] in ['CF'] ) and ( parts[5] in ['gid'] ):
                out['SNP_CN_4407927_T276G_E92D_gid'] = "1"

    for variant in variants:
        output2.write("," + str(int(variant in out)))
    output2.write("\n")

if __name__ == '__main__':
    main()
