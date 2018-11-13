#!/usr/bin/env python3
#Octobre 2018
#Gravey François

import subprocess
import os
import re
from argparse import ArgumentParser
from shutil import copyfile


input_file = '/Users/Francois/Documents/projets/ecloacae/patric_db/analyses/\
blast/romr_roma/clustering/romR_clustering.fasta.clstr'
gene = 'romR'
outputdir = '/Users/Francois/Documents/projets/ecloacae/patric_db/analyses/\
blast/romr_roma/'
filename = 'romr_clustering'

regex = re.compile("^>")
dico = {}
tout = []
with open(input_file, 'r') as filin:
    for ligne in filin:
        if regex.search(ligne):
            dico[ligne[:-1]] = []

with open(input_file, 'r') as filin:
    for ligne in filin:
        tout.append(ligne[:-1])

for cluster in dico.keys():
    i = int(tout.index(cluster)) +1
    while regex.search(tout[i]) == None:
        if (int(i) < int((len(tout)-1))) == True:
            dico[cluster].append(tout[i].split('\t')[1].split(",")[1].split("_")[0].replace(">",""))
            i = i+1
        else:
            dico[cluster].append(tout[i].split('\t')[1].split(",")[1].split("_")[0].replace(">",""))
            break

with open()
for cluster in dico.keys():
    for fasta in dico[cluster]:
        print("{};{}{}".format(fasta, gene,cluster.replace(">", "_")))
#print(int(i))
##print(len(tout)-1)
#print(
