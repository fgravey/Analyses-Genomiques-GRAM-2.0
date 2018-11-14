#!/usr/bin/env python3
## Octobre 2018
### Gravey François

import re

class ReturnValue(object):
  def __init__(self, position, information):
     self.position = position
     self.information = information



#Regex in order to not include CDS into the research
regex_3 = re.compile('CDS')

def read_vcf(file):
    ## Regex definition
    #regex specific to header lines in the beguining og the vcf file
    regex = re.compile('^"')
    regex_2 = re.compile('^#')

    #Container definition
    position = []
    information = {}

    with open(file, 'r') as filin:
        for ligne in filin:
            #Parsing the file and do not take all the headers lines in the beguining of the file
            if regex.search(ligne) == None and regex_2.search(ligne) == None:
                # Variable definition
                pos = ligne.split('\t')[1] #position in the genome of the ref
                nt_ref = ligne.split('\t')[3] #nucleotide in the ref genome
                nt_strain = ligne.split('\t')[4] #nucleotide in the genome of the strain
                qual = ligne.split('\t')[5] #Quality score of the mapping
                info = ligne.split('\t')[7] #information in the .vcf file about the SNP calling
                form = ligne.split('\t')[8] #format header in the .vcf file

                position.append("{}:{}:{}".format(pos,nt_ref,nt_strain))

                information[pos] = "{},{},{},{},{}".format(nt_ref,nt_strain,qual,info,form)

    #End of function
    return ReturnValue(position,information)

def discordance(file1,file2):
    """Look for snp with are present into the file1 but not into the file2
    regarding the position of the nucleotides and the nature of the substitutions"""

    # Container definition
    discor = [] #list which will contain all the snp fond in the file one but not in the file2

    for pos in read_vcf(file1).position:
        if pos in read_vcf(file2).position:
            continue
    else:
        discor.append("{}".format(pos))

    #End of the function
    return(discor)

def verif_snp(file1,file2):
    #Variable definition
    liste1 = discordance(file1,file2) #Contained all the snp in the file1 but not in the file2
    liste2 = discordance(file2,file1) #Contained all the snp in the file2 but not in the file1

    if len(liste1) == 0 and len(liste2) == 0:
        for snp in liste1:
            if snp in liste2:
                print("ALERTE FILRTE FAILED")
                break

        for snp in liste2:
            if snp in liste1:
                print("ALERTE FILRTE FAILED")
                break

# distance = (len(discordance_1) + len(discordance_2))

# bed = []
# with open("/Users/Francois/Desktop/CP012165_1.bed", "r") as filin:
    # for ligne in filin:
        # if regex_3.search(ligne) == None:
            #print(ligne.split("\t"))
            # info = [ligne.split("\t")[1], ligne.split("\t")[2], ligne.split("\t")[7], ligne.split("\t")[9]]
            # bed.append(info)
#
#
# snp_interets_1 = []
# for snp in discordance_1:
    # for pos in bed:
        # if int(pos[0]) <= int(snp.split(":")[0]) and int(pos[1]) > int(snp.split(":")[0]):
            #variable definition
            # snp_position = snp.split(":")[0]
            # snp_information = dico_1[snp.split(":")[0]]
            # gene_start = pos[0]
            # gene_end = pos[1]
            # type = pos[2]
            # gene_info = pos[3]
#
            #Adding all the snp of interest into the snp_interets_1 list
            # snp_interets_1.append("SNP,{},{},{},{},{},{}".format(snp_position,\
            # snp_information,gene_start, gene_end, type,gene_info))
#
file1 = "/Users/Francois/Desktop/FAY1_test.vcf"
file2 = "/Users/Francois/Desktop/FAY2_test.vcf"

print(discordance(file1,file2))
print(verif_snp(file1,file2))

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
