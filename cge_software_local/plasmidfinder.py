#!/usr/bin/env python3
## Septembre 2018
### Gravey François

### module Loading
import subprocess
import re
from argparse import ArgumentParser

##########################################################################
# PARSE COMMAND LINE OPTIONS
##########################################################################

parser = ArgumentParser()
parser.add_argument("-l", "--list", dest="list",\
help="list which contains all the names of the strains", default='')
parser.add_argument("-f", "--fasta_dir", dest="fasta_dir",help="", \
default='Directory which contains all the fasta files')
parser.add_argument("-e", "--extension", dest="extension",\
help="Any informations such as .agp.fasta .awked.fa .scfd.fasta .fasta .fa", default='')
parser.add_argument("-o", "--outputPath", dest="out_path",help="Path to blast output", default='')
args = parser.parse_args()

# Defining varibales

liste = args.list
fasta_dir = args.fasta_dir
outputdir = args.out_path
fasta_extension = args.extension

#### unchanging path
plasmidfinder = '/Users/Francois/cge_softwares/plasmidfinder/plasmidfinder.py'
database = '/Users/Francois/cge_data_bases/plasmidfinder_db'

#listing all the files which will be working on
travail = []
with open(liste, 'r') as filin:
    for nom in filin:
        travail.append(nom[:-1])

print("Voici la liste des fichiers sur lesquels vous allez travailler", travail)
print("Le nombre de fichier est de : {}".format(len(travail)))

resultats = ["Strain;Database;Plasmide;Identities;Alignment.Length;Template.Length;\
Position.in.reference;Contig;Position.in.contig;Note;Accession no.\n"]
for nom in travail:
    fasta = "{}/{}{}".format(fasta_dir, nom, fasta_extension)
    outputdir2 = "{}/{}".format(outputdir,nom)
    subprocess.run(["mkdir", outputdir2])
    print("##################################################")
    print("working on {}.{} file".format(nom, fasta_extension))
    print("##################################################")
    subprocess.run(["python", "{}".format(plasmidfinder), "-p", "{}".format(database),\
     "-i", "{}".format(fasta), "-o", "{}".format(outputdir2, "-d", "Enterobacteriaceae")])
    plasmide = []
    with open("{}/results_tab.txt".format(outputdir2),"r") as filin:
        for ligne in filin:
            plasmide.append(ligne[:-1])
        if len(plasmide) == 1:
            resultats.append("{};{};no plasmid found;-;-;-;-;-;-;-;-\n"\
            .format(nom,"Enterobacteriaceae"))
        else:
            for i in range(1,len(plasmide)):
                resultats.append("{};{}\n".format(nom,plasmide[i].replace("\t",";")))

with open("{}/plasmidfinder_results.csv".format(outputdir), "w") as filout:
    for res in resultats:
        filout.write(res)
