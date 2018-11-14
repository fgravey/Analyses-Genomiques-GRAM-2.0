#!/usr/bin/env python3
## Octobre 2018
### Gravey François

import re

## Regex definition
#regex specific to header lines in the beguining og the vcf file
regex = re.compile('^"')
#regex specific to header lines in the beguining og the vcf file
regex_2 = re.compile('^#')
#Regex in order to not include CDS into the research
regex_3 = re.compile('CDS')

#Container definition
pos_vcf_1 = []
pos_vcf_2 = []
dico_1 = {}
dico_2 = {}
discordance_1 = []
discordance_2 = []

with open("/Users/Francois/Desktop/FAY1_test.vcf", 'r') as filin:
    for ligne in filin:
        #Parsing the file and do not take all the headers lines in the beguining of the file
        if regex.search(ligne) == None and regex_2.search(ligne) == None:
            # for each position in the vcf file adding it in the pos_vcf_1 list
            pos_vcf_1.append("{}:{}:{}".format(ligne.split('\t')[1],\
            ligne.split('\t')[3],ligne.split('\t')[4]))
            dico_1[ligne.split('\t')[1]] = "{},{},{},{},{}".format(\
            ligne.split('\t')[3],ligne.split('\t')[4],ligne.split('\t')[5],\
            ligne.split('\t')[7],ligne.split('\t')[8])
#

with open("/Users/Francois/Desktop/FAY2_test.vcf", 'r') as filin:
    for ligne in filin:
        #Parsing the file and do not take all the headers lines in the beguining of the file
        if regex.search(ligne) == None and regex_2.search(ligne) == None:
            # for each position in the vcf file adding it in the pos_vcf_1 list
            pos_vcf_2.append("{}:{}:{}".format(ligne.split('\t')[1],\
            ligne.split('\t')[3],ligne.split('\t')[4]))
#
print(dico_1)
print(pos_vcf_1)
print(pos_vcf_2)

for pos in pos_vcf_1:
    if pos in pos_vcf_2:
        continue
    else:
        discordance_1.append("{}".format(pos))

for pos in pos_vcf_2:
    if pos in pos_vcf_1:
        continue
    else:
        discordance_2.append("{}".format(pos))
#
print(discordance_1)
print(discordance_2)
distance = (len(discordance_1) + len(discordance_2))
print(distance)
# print(len(discordance_1))
# print(len(discordance_2))
#
# bed = []
# with open("/Users/Francois/Desktop/CP012165_1.bed", "r") as filin:
    # for ligne in filin:
        # if regex_3.search(ligne) == None:
            #print(ligne.split("\t"))
            # info = [ligne.split("\t")[1], ligne.split("\t")[2], ligne.split("\t")[7], ligne.split("\t")[9]]
            # bed.append(info)
#
# genes_interets_1 = []
# for snp in discordance_1:
    # for pos in bed:
        # if int(pos[0]) <= int(snp) and int(pos[1]) > int(snp):
            # genes_interets_1.append(pos)
#
#
# genes_interets_2 = []
# for snp in discordance_2:
    # for pos in bed:
        # if int(pos[0]) <= int(snp) and int(pos[1]) > int(snp):
            # genes_interets_2.append(pos)
#
# resultats_1 = []
# for gene in genes_interets_1:
    # if gene in resultats_1:
        # continue
    # else:
        # resultats_1.append(gene)
#
#
# resultats_2 = []
# for gene in genes_interets_2:
    # if gene in resultats_2:
        # continue
    # else:
        # resultats_2.append(gene)
#
#
# with open("/Users/Francois/Desktop/SNP_annotation_FAY1.csv", "w") as filout:
    # for ligne in resultats_1:
        # filout.write(";".join(ligne))
#
# with open("/Users/Francois/Desktop/SNP_annotation_FAY9.csv", "w") as filout:
    # for ligne in resultats_2:
        # filout.write(";".join(ligne))
