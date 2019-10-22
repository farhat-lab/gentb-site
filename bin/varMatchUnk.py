#!/usr/bin/env python

import argparse
import json
import re
import sys

parser = argparse.ArgumentParser(
    description="processes prediction json and var and updates json with other "
    "mutations to display in the prediction report, and includes a "
    "new category of lineage mutations. Includes a set of drugs and"
    "gene names"
)

parser.add_argument("var_file", help="var file")
parser.add_argument(
    "lineage_snps_file",
    help="file with snps that are lineage specific and not resistance causing from WALKER paper",
)
parser.add_argument(
    "json_file",
    help="file that contains R prediction and list of important and other mutations",
)

args = parser.parse_args()
name = args.var_file

drugs = set()

drug_mapping = {
    "AMK": "AMIKACIN",
    "CAP": "CAPREOMYCIN",
    "CIP": "CIPROFLOXACIN",
    "EMB": "ETHAMBUTOL",
    "ETH": "ETHIONAMIDE",
    "INH": "ISONIAZID",
    "KAN": "KANAMYCIN",
    "LEVO": "LEVOFLOXACIN",
    "OFLX": "OFLOXACIN",
    "PAS": "PARA-AMINOSALICYLIC_ACID",
    "PZA": "PYRAZINAMIDE",
    "RIF": "RIFAMPICIN",
    "STR": "STREPTOMYCIN",
}

drugs_to_genes = {
    "AMK": set(["rrs"]),
    "CAP": set(["rrs", "tlyA"]),
    "CIP": set(["gyrA", "gyrB"]),
    "EMB": set(["embB", "embC", "embA", "Rv3806c", "promoter-embA-embB"]),
    "ETH": set(["ethA", "promoter-fabG1-inhA", "inhA"]),
    "INH": set(
        [
            "katG",
            "promoter-fabG1-inhA",
            "inhA",
            "kasA",
            "promoter-ahpC",
            "promoter-embA-embB",
        ]
    ),
    "KAN": set(["rrs", "rrl", "eis", "inter-eis-Rv2417c"]),
    "LEVO": set(["gyrA", "gyrB"]),
    "OFLX": set(["gyrA", "gyrB"]),
    "PAS": set(["thyA", "inter-thyX-hsdS.1", "folC"]),
    "PZA": set(["pncA", "promoter-pncA", "rpsA"]),
    "RIF": set(["rpoB"]),
    "STR": set(["rpsL", "gid", "rrs", "inter-rrs-rrl"]),
}


class Variant(object):
    def __init__(self, gene_name, codon_location, AA_change,
                 name=None, test_name=None, drug=None):
        self.gene_name = gene_name
        self.codon_location = codon_location
        self.AA_change = AA_change
        self.name = name
        self.drug = drug
        # print(test_name)

    def compare_variant(self, variant):
        return (self.gene_name == variant.gene_name) and (
            self.codon_location == variant.codon_location
        )  ##and (self.AA_change == variant.AA_change) ##will include the AA matching below for coding snps and Indel intergenic/promoters

    def __str__(self):
        if self.AA_change:
            return_string = self.AA_change[0] + self.codon_location + self.AA_change[1]
        else:
            return_string = self.codon_location
        return "{}_{}".format(self.gene_name, return_string)


"""
Example piece of data in json file, major discrepancy with current var are promoter and intergenic regions and indels:

SNP_CN_2155168_C944G_S315T_katG
SNP_P_1673425_C15T_promoter_fabG1.inhA
DEL_I_1476813_d86GGGAG_inter_rrl_rrf

for testing

imp['INH'] = ['SNP_CN_2155168_C944G_S315T_katG']
imp['EMB'] = ['SNP_CN_4247431_G918C_M306I_embB']
imp['PAS'] = ['SNP_CN_3073868_T604C_T202A_thyA','INS_P_3074519_i48G_promoter_thyA']
imp['STR'] = ['SNP_CZ_4407731_G472A_R158._gid','INS_I_1471827_i19AGA_inter_murA_rrs']
imp['ETH'] = ['DEL_CD_4326366_d1108TGTAGGCCATCG_370_ethA']
imp['PZA'] = ['SNP_P_2289253_A12C_promoter_pncA']
imp['RIF'] = ['INS_CI_761103_i1296TTC_433_rpoB']
imp['LEVO'] = ['SNP_I_7268_C34T_inter_gyrB_gyrA','SNP_P_5075_C48T_promoter_gyrB.gyrA','INS_P_5079_i44G_promoter_gyrB.gyrA', 'INS_CF_9373_i2071T_691_gyrA']

"""
with open(args.json_file, "r") as f:
    datastore = json.load(f)

drs = [
    "INH",
    "RIF",
    "PZA",
    "EMB",
    "STR",
    "ETH",
    "KAN",
    "CAP",
    "AMK",
    "CIP",
    "LEVO",
    "OFLX",
    "PAS",
]

imp = {}
oth = {}

#import sys

j = 0
for d in drs:
    imp[d] = []
    oth[d] = []
    for i in range(0, 5):
        imp[d].append(datastore[1][datastore[1].keys()[0]][i][j])
        oth[d].append(datastore[2][datastore[2].keys()[0]][i][j])
    j = j + 1


imp_variants_identified = []
oth_variants_identified = []

"""

can use the code below to process all randomforest variants variant_name_list.csv if needed

"""

for d in drs:
    for mut in imp[d]:
        if mut:
            #sys.stderr.write(mut)
            #sys.stderr.write('\n')
            type_change_info = mut.split("_")
            if type_change_info[0] == "SNP":
                if type_change_info[1] == "CN":
                    gene_name, codonAA = type_change_info[5], type_change_info[4]
                    type_change, codon_position = (
                        codonAA[0] + codonAA[len(codonAA) - 1],
                        codonAA[1 : len(codonAA) - 1],
                    )
                elif type_change_info[1] == "CZ":
                    gene_name, codonAA = (
                        type_change_info[5],
                        type_change_info[4].replace(".", "*"),
                    )
                    type_change, codon_position = (
                        codonAA[0] + codonAA[len(codonAA) - 1],
                        codonAA[1 : len(codonAA) - 1],
                    )
                elif type_change_info[1] == "P":
                    gene_name, codonAA = type_change_info[5], type_change_info[3]
                    if "." in gene_name:
                        gene_name = "promoter-" + gene_name.replace(".", "-")
                    else:
                        gene_name = "promoter-" + gene_name
                    type_change, codon_position = (
                        codonAA[0] + codonAA[len(codonAA) - 1],
                        codonAA[1 : len(codonAA) - 1],
                    )
                elif type_change_info[1] == "N":
                    gene_name, codonAA = type_change_info[4], type_change_info[3]
                    type_change, codon_position = (
                        codonAA[0] + codonAA[len(codonAA) - 1],
                        codonAA[1 : len(codonAA) - 1],
                    )
                elif type_change_info[1] == "I":
                    gene_name = (
                        "inter-" + type_change_info[5] + "-" + type_change_info[6]
                    )
                    codonAA = type_change_info[3]
                    type_change, codon_position = (
                        codonAA[0] + codonAA[len(codonAA) - 1],
                        codonAA[1 : len(codonAA) - 1],
                    )
                    if type_change_info[5] == "gyrB":
                        gene_name = "promoter-gyrB-gyrA"
            elif type_change_info[0] in ["INS", "DEL"] and type_change_info[1] == "CF":
                gene_name, codon_position = type_change_info[5], type_change_info[4]
                codonAA = type_change_info[3]
                m = re.search(r"[ACGT]+", codonAA)
                if m:
                    indel_seq = m.group()
                type_change = (
                    codonAA[0] + type_change_info[1][1] + indel_seq.replace("\n", "")
                )  # e.g dFG or iFGA
            elif type_change_info[0] in ["INS", "DEL"] and type_change_info[1] in [
                "CD",
                "CI",
            ]:
                gene_name, codon_position = type_change_info[5], type_change_info[4]
                codonAA = type_change_info[3]
                m = re.search(r"[ACGT]+", codonAA)
                if m:
                    indel_seq = m.group()
                type_change = (
                    codonAA[0] + type_change_info[1][1] + indel_seq.replace("\n", "")
                )  # e.g dIG or iDGA
            elif type_change_info[0] in [
                "INS",
                "DEL",
            ]:  # must be indel in promoter or intergenic region
                if type_change_info[1] == "P":
                    gene_name, codonAA = type_change_info[5], type_change_info[3]
                    if "." in gene_name:
                        gene_name = "promoter-" + gene_name.replace(".", "-")
                    else:
                        gene_name = "promoter-" + gene_name
                elif type_change_info[1] == "I":
                    gene_name = (
                        "inter-" + type_change_info[5] + "-" + type_change_info[6]
                    )
                    codonAA = type_change_info[3]
                    m = re.search(r"(\d+)([ACGT]+)", codonAA)
                    if m:
                        codon_position = m.group(1)
                        indel_seq = m.group(2)
                    type_change = codonAA[0] + indel_seq.replace("\n", "")
            test = Variant(gene_name, codon_position, type_change)
            imp_variants_identified.append(test)
            drugs_to_genes[d].add(gene_name)

for d in drs:
    for mut in oth[d]:
        if mut:
            #sys.stderr.write(mut)
            #sys.stderr.write('\n')
            type_change_info = mut.split("_")
            if type_change_info[0] == "SNP":
                if type_change_info[1] == "CN":
                    gene_name, codonAA = type_change_info[5], type_change_info[4]
                    type_change, codon_position = (
                        codonAA[0] + codonAA[len(codonAA) - 1],
                        codonAA[1 : len(codonAA) - 1],
                    )
                elif type_change_info[1] == "CZ":
                    gene_name, codonAA = (
                        type_change_info[5],
                        type_change_info[4].replace(".", "*"),
                    )
                    type_change, codon_position = (
                        codonAA[0] + codonAA[len(codonAA) - 1],
                        codonAA[1 : len(codonAA) - 1],
                    )
                elif type_change_info[1] == "P":
                    gene_name, codonAA = type_change_info[5], type_change_info[3]
                    if "." in gene_name:
                        gene_name = "promoter-" + gene_name.replace(".", "-")
                    else:
                        gene_name = "promoter-" + gene_name
                    type_change, codon_position = (
                       	codonAA[0] + codonAA[len(codonAA) - 1],
                        codonAA[1 : len(codonAA) - 1],
                    )
                elif type_change_info[1] == "N":
                    gene_name, codonAA = type_change_info[4], type_change_info[3]
                    type_change, codon_position = (
                       	codonAA[0] + codonAA[len(codonAA) - 1],
                        codonAA[1 : len(codonAA) - 1],
                    )
                elif type_change_info[1] == "I":
                    gene_name = (
                       	"inter-" + type_change_info[5] + "-" + type_change_info[6]
                    )
                    codonAA = type_change_info[3]
                    type_change, codon_position = (
                       	codonAA[0] + codonAA[len(codonAA) - 1],
                        codonAA[1 : len(codonAA) - 1],
                    )
                    if type_change_info[5] == "gyrB":
                       	gene_name = "promoter-gyrB-gyrA"
            elif type_change_info[0] in ["INS", "DEL"] and type_change_info[1] == "CF":
                gene_name, codon_position = type_change_info[5], type_change_info[4]
                codonAA = type_change_info[3]
                m = re.search(r"[ACGT]+", codonAA)
                if m:
                    indel_seq = m.group()
                type_change = (
                    codonAA[0] + type_change_info[1][1] + indel_seq.replace("\n", "")
               	)  # e.g dFG or iFGA
            elif type_change_info[0] in ["INS", "DEL"] and type_change_info[1] in [
                "CD",
                "CI",
            ]:
                gene_name, codon_position = type_change_info[5], type_change_info[4]
                codonAA = type_change_info[3]
                m = re.search(r"[ACGT]+", codonAA)
                if m:
                    indel_seq = m.group()
                type_change = (
                    codonAA[0] + type_change_info[1][1] + indel_seq.replace("\n", "")
                )  # e.g dIG or iDGA
            elif type_change_info[0] in [
               	"INS",
                "DEL",
            ]:  # must be indel in promoter or intergenic region
                if type_change_info[1] == "P":
                    gene_name, codonAA = type_change_info[5], type_change_info[3]
                    if "." in gene_name:
                        gene_name = "promoter-" + gene_name.replace(".", "-")
                    else:
                        gene_name = "promoter-" + gene_name
                elif type_change_info[1] == "I":
                    gene_name = (
                        "inter-" + type_change_info[5] + "-" + type_change_info[6]
                    )
                    codonAA = type_change_info[3]
                    type_change, codon_position = (
                       	codonAA[0] + codonAA[len(codonAA) - 1],
                       	codonAA[1 : len(codonAA) - 1],
                    )
                    if type_change_info[5] == "gyrB":
                        gene_name = "promoter-gyrB-gyrA"
            elif type_change_info[0] in ["INS", "DEL"] and type_change_info[1] == "CF":
               	gene_name, codon_position = type_change_info[5], type_change_info[4]
                codonAA = type_change_info[3]
                m = re.search(r"[ACGT]+", codonAA)
                if m:
                    indel_seq = m.group()
                type_change = (
                    codonAA[0] + type_change_info[1][1] + indel_seq.replace("\n", "")
                )  # e.g dFG or iFGA
            elif type_change_info[0] in ["INS", "DEL"] and type_change_info[1] in [
                "CD",
                "CI",
            ]:
              	gene_name, codon_position = type_change_info[5], type_change_info[4]
                codonAA = type_change_info[3]
                m = re.search(r"[ACGT]+", codonAA)
                if m:
                    indel_seq = m.group()
                type_change = (
                    codonAA[0] + type_change_info[1][1] + indel_seq.replace("\n", "")
                )  # e.g dIG or iDGA
            elif type_change_info[0] in [
                "INS",
                "DEL",
            ]:  # must be indel in promoter or intergenic region
                if type_change_info[1] == "P":
                    gene_name, codonAA = type_change_info[5], type_change_info[3]
                    if "." in gene_name:
                        gene_name = "promoter-" + gene_name.replace(".", "-")
                    else:
                        gene_name = "promoter-" + gene_name
                elif type_change_info[1] == "I":
                    gene_name = (
                        "inter-" + type_change_info[5] + "-" + type_change_info[6]
                    )
                    codonAA = type_change_info[3]
                    m = re.search(r"(\d+)([ACGT]+)", codonAA)
                    if m:
                        codon_position = m.group(1)
                       	indel_seq = m.group(2)
                    type_change = codonAA[0] + indel_seq.replace("\n", "")
            test = Variant(gene_name, codon_position, type_change)
            oth_variants_identified.append(test)
            drugs_to_genes[d].add(gene_name)



# sys.stderr.write("\t".join(drugs_to_genes['INH']))

"""
Now break down lineage specific snps, in same format as var, but genomic coords are placeholder

"""

lineage_snps = []

for line in open(args.lineage_snps_file, "r").readlines()[1:]:
    line = line.replace("\n", "")
    type_change_info = line.split("_")
    # sys.stderr.write(line+"\n")
    if type_change_info[1] == "CN":  # not currently used or tested
        gene_name, codonAA = type_change_info[5], type_change_info[4]
        type_change, codon_position = (
            codonAA[0] + codonAA[len(codonAA) - 1],
            codonAA[1 : len(codonAA) - 1],
        )
    elif type_change_info[1] == "P":
        gene_name, codonAA = type_change_info[4], type_change_info[3]
        type_change, codon_position = (
            codonAA[0] + codonAA[len(codonAA) - 1],
            codonAA[1 : len(codonAA) - 1],
        )
    elif type_change_info[1] == "N":
        gene_name, codonAA = type_change_info[4], type_change_info[3]
        type_change, codon_position = (
            codonAA[0] + codonAA[len(codonAA) - 1],
            codonAA[1 : len(codonAA) - 1],
        )
    elif type_change_info[1] == "I":  # not currently used or tested
        gene_name, codonAA = type_change_info[4], type_change_info[3]
        type_change, codon_position = (
            codonAA[0] + codonAA[len(codonAA) - 1],
            codonAA[1 : len(codonAA) - 1],
        )
    elif (
        type_change_info[0] == "INS" and type_change_info == "F"
    ):  # not currently used or tested
        gene_name, codon_position = type_change_info[4], type_change_info[3]
        type_change, codon_position = (
            codonAA[0] + codonAA[len(codonAA) - 1],
            codonAA[1 : len(codonAA) - 1].replace("\n", ""),
        )
    elif (
        type_change_info[0] == "DEL" and type_change_info[1] == "F"
    ):  # not currently used or tested
        gene_name, codon_position = type_change_info[4], type_change_info[3]
        type_change, codon_position = (
            codonAA[0] + codonAA[len(codonAA) - 1],
            codonAA[1 : len(codonAA) - 1].replace("\n", ""),
        )

    test = Variant(gene_name, codon_position, type_change)
    lineage_snps.append(test)


""" 
SNP_CN_7585_G284C_S95T_gyrA

SNP_P_781395_T165C_promoter-rpsL

SNP_N_1474571_G914A_rrl

INS_CF_2289050_i192A_65I_pncA

"""

all = []
with open(args.var_file, "r") as f:
    for line in f:
        # sys.stderr.write(line)
        result = line.rstrip().split("\t")
        all.append(result[5])
        # sys.stderr.write(result[5]+"\n")

new_variants_identified = {}
lineage_variants_identified = {}
moth={}

for d in drs:
    new_variants_identified[d] = []
    lineage_variants_identified[d] = []
    moth[d] = []
    # sys.stderr.write(d+':\t')
    for i in drugs_to_genes[d]:
        # sys.stderr.write(i+'\t')
        pattern = i
        pattern = r"\s?([\w_\.]+%s[\w_\.]*)\s?" % pattern
        muts = tuple(re.finditer(pattern, ",".join(all))) #, re.IGNORECASE)) #respect lowercase
        if len(muts) >= 1:
            #sys.stderr.write("with pattern "+i+" found "+str(len(muts))+" variant(s) for drug "+d+"\n")
            for j in range(0, len(muts)):
                var = muts[j].group()
                #sys.stderr.write(str(var)+"\n")
                type_change_info = var.split("_")
                sys.stderr.write("variant is "+type_change_info[0]+"\n")
                if type_change_info[1] != "CS":
                    if type_change_info[0] == "SNP" and type_change_info[1] in [
                        "CN",
                        "CZ",
                    ]:
                        gene_name, codonAA = type_change_info[5], type_change_info[4]
                        type_change, codon_position = (
                            codonAA[0] + codonAA[len(codonAA) - 1],
                            codonAA[1 : len(codonAA) - 1],
                        )
                    elif type_change_info[0] == "SNP" and type_change_info[1] == "P":
                        gene_name, codonAA = type_change_info[4], type_change_info[3]
                        type_change, codon_position = (
                            codonAA[0] + codonAA[len(codonAA) - 1],
                            codonAA[1 : len(codonAA) - 1],
                        )
                    elif type_change_info[0] == "SNP" and type_change_info[1] == "N":
                        gene_name, codonAA = type_change_info[4], type_change_info[3]
                        type_change, codon_position = (
                            codonAA[0] + codonAA[len(codonAA) - 1],
                            codonAA[1 : len(codonAA) - 1],
                        )
                    elif type_change_info[0] == "SNP" and type_change_info[1] == "I":
                        gene_name, codonAA = type_change_info[4], type_change_info[3]
                        type_change, codon_position = (
                            codonAA[0] + codonAA[len(codonAA) - 1],
                            codonAA[1 : len(codonAA) - 1],
                        )
                    elif type_change_info[0] in ["INS", "DEL"] and type_change_info[
                        1
                    ] in ["CF", "CD", "CI"]:
                        gene_name, codon_position = (
                            type_change_info[5],
                            type_change_info[4],
                        )
                        if type_change_info[0] == "INS":
                            m = re.search(r"\d+", codon_position)
                            if m:
                                codon_position = m.group()
                        codonAA = type_change_info[3]
                        m = re.search(r"[ACGT]+", codonAA)
                        if m:
                            indel_seq = m.group()
                        type_change = (
                            codonAA[0]
                            + type_change_info[1][1]
                            + indel_seq.replace("\n", "")
                        )  # e.g dFG or iFGA
                    elif type_change_info[0] in [
                        "INS",
                        "DEL",
                    ]:  # must be indel in promoter or intergenic region
                        gene_name, codonAA = type_change_info[4], type_change_info[3]
                        codonAA = type_change_info[3]
                        m = re.search(r"(\d+)([ACGT]+)", codonAA)
                        if m:
                            codon_position = m.group(1)
                            indel_seq = m.group(2)
                        type_change = codonAA[0] + indel_seq.replace("\n", "")
                    test = Variant(
                        gene_name, codon_position, type_change, name=str(var)
                    )
                    lineage_test = 0
                    oth_test=0
                    imp_test=0
                    for lin in lineage_snps:
                        if lin.compare_variant(test): #codon_location and gene_names are the same
                            if lin.AA_change == test.AA_change:
                                lineage_test = 1
                    for imp in imp_variants_identified:
                        if imp.compare_variant(test):
                            if imp.AA_change == test.AA_change:
                                imp_test = 1
                            elif imp.AA_change[0:2] in ["dF","iF"]:
                                imp_test = 1 
                    for oth in oth_variants_identified:
                        if oth.compare_variant(test):
                            if oth.AA_change == test.AA_change:
                                oth_test = 1
                            elif oth.AA_change[0:2] in ["dF", "iF"]:
                                oth_test = 1
                    if oth_test == 1 and lineage_test == 1:
                        lineage_variants_identified[d].append(var)
                    if oth_test == 1 and lineage_test == 0:
                        #sys.stderr.write("other variant found and is "+ var+'\n')
                        moth[d].append(str(var))
                    if imp_test == 0 and oth_test == 0 and lineage_test == 0:
                        new_variants_identified[d].append(var)


results1 = {}
results2 = {}
other=datastore[2]
v1 = []
v2 = []
for d in drs:
    v1.append(len(new_variants_identified[d]))
    v2.append(len(lineage_variants_identified[d]))
l1 = max(v1)
l2 = max(v2)

if l1 == 0:
    results1[0] = ["Null"] * len(drs)
else:
    for i in range(0, l1):
        results1[i] = []
        results1[i] = ["Null"] * len(drs)
        j = -1
	# sys.stderr.write(d+' out\n')
        for d in range(0, len(drs)):
            j = j + 1
            if v1[d] > i:
                if len(new_variants_identified[drs[d]]) > 0:
                    # sys.stderr.write(d+' in!\n')
                    results1[i][j] = new_variants_identified[drs[d]][i]

if l2 == 0:
    results2[0] = ["Null"] * len(drs)
else:
    for i in range(0, l2):
        results2[i] = []
        results2[i] = ["Null"] * len(drs)
        j = -1
	# sys.stderr.write(d+' out\n')
        for d in range(0, len(drs)):
            j = j + 1
            if v2[d] > i:
                if len(lineage_variants_identified[drs[d]]) > 0:
                    # sys.stderr.write(d+' in!\n')
                    results2[i][j] = lineage_variants_identified[drs[d]][i]

for i in range(0, 5):
        j = -1
        for d in range(0, len(drs)):
            j = j + 1
            if len(moth[drs[d]][i:]) > 0:
                    other.items()[0][1][i][j] = moth[drs[d]][i]
            else:
                    other.items()[0][1][i][j] = None

def append_results(filename, *results):
    """Write the resulting structures to the output filename"""
    with open(filename, "r") as infile:
        structure = json.loads(infile.read())
        structure.pop()
        structure.extend(results)

    with open(filename, "w") as outfile:
        outfile.write(json.dumps(structure, indent=2))


append_results(args.json_file, other, results1, results2)
