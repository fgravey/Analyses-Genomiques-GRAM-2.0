#!/usr/bin/env python3
## Septembre 2018
### Gravey François

### module Loading
import subprocess
import re
import glob
from argparse import ArgumentParser


# PARSE COMMAND LINE OPTIONS
##########################################################################
parser = ArgumentParser()
parser.add_argument("-l", "--list", dest="list", \
help="list which contains all the name of the strains", default='')
parser.add_argument("-o", "--outputPath", dest="out_path",\
help="Path to blast output", default='')
parser.add_argument("-f", "--fasta_path", dest="fasta_path_dir",\
help="Path to the directory which contains all the fasta to analyse", default='')
parser.add_argument("-e", "--extension", dest="extension",\
help="fasta extension such as .fasta .fa .awked.fasta .agp.fasta ", default='')
args = parser.parse_args()

###############################################################################
# Variables difinition
liste = args.list
fasta_dir = args.fasta_path_dir
outputdir = args.out_path
fasta_extension = args.extension

###############################################################################
#listing all the files which will be working on
travail = []
with open(liste, 'r') as filin:
    for nom in filin:
        travail.append(nom[:-1])

print("Voici la liste des fichiers sur lesquels vous allez travailler", travail)
print("Le nombre de fichier est de : {}".format(len(travail)))

###############################################################################
#launching rgi software
for nom in travail:
    fasta = "{}/{}{}".format(fasta_dir, nom, fasta_extension)
    subprocess.run(["rgi", "main", "--input_sequence", "{}".format(fasta), "--output_file",\
     "{}/{}".format(outputdir,nom), "--input_type", "contig"])

###############################################################################
# Cleaning the temporary files
for fichier in glob.glob("{}/*.temp*".format(outputdir)):
    subprocess.run(["rm", fichier])

###############################################################################
# Creating the heatmap and clustering the strains based on atb resistance genes
subprocess.run(["rgi", "heatmap", "-i", "{}".format(outputdir), "-clus", "samples"])
