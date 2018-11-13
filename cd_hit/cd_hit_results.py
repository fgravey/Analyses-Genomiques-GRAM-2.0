#!/usr/bin/env python3
#Octobre 2018
#Gravey François

# Loading modules
import os
import re
from argparse import ArgumentParser

def parse_cd_hit(input_file):
    #Regex definition
    regex = re.compile("^>")

    #Containers definition
    dico = {} ## Which will contain the name of all the clusters recevered
    resultats = [] ## Which will contain the information in the input_file

    #Reading the input_file and each cluster became a key of the dico dictionnary
    with open(input_file, 'r') as filin:
        for ligne in filin:
            if regex.search(ligne):
                dico[ligne[:-1]] = []

    #Reading the input_file and all the information available is put on the resultats container
    with open(input_file, 'r') as filin:
        for ligne in filin:
            resultats.append(ligne[:-1])

    # For each key in dico == for each cluster find by cd-hit looking for all the fasta
    for cluster in dico.keys():
        i = int(resultats.index(cluster)) +1
        while regex.search(resultats[i]) == None:
            # Mandatory in order to no be out of range
            # For each key, all the fasta are added into the dico[key] dictionnary
            if (int(i) < int((len(resultats)-1))) == True:
                dico[cluster].append(resultats[i].split('\t')[1].split(",")[1].split("_")[0].replace(">",""))
                i = i+1
            else:
                dico[cluster].append(resultats[i].split('\t')[1].split(",")[1].split("_")[0].replace(">",""))
                break
    #End of function
    return(dico)

if __name__ == "__main__":

    # Parse command line options
    ##########################################################################
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_file", dest="input_file",\
    help="absolute path to the cd hits .clstr file", default='')
    parser.add_argument("-g", "--gene", dest="gene",\
    help="gene you are working on", default='')
    parser.add_argument("-o", "--outputPath", dest="out_path",\
    help="Path to putt the .csv file", default='')
    parser.add_argument("-filename", "--filename", dest="filename",\
    help="summary output file name", default='clustering_cd_hits')
    args = parser.parse_args()

    ###############################################################################
    # Variables difinition
    input_file = args.input_file
    gene = args.gene
    outputdir = args.out_path
    filename = args.filename

    with open("{}{}.csv".format(outputdir, filename), 'w') as filout:
        filout.write("Souche;{}_cd-hits_cluster\n".format(gene))
        for cluster in parse_cd_hit(input_file).keys():
            for fasta in parse_cd_hit(input_file)[cluster]:
                filout.write("{};{}_{}\n".format(fasta[1:].replace(" ",""),\
                 gene,cluster.replace(">", "").replace(" ","_")))
