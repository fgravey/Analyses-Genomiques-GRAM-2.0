#!/usr/bin/env python3
## Septembre 2018
### Gravey François

### module Loading
import subprocess
import re

##### changing path
liste = '/Volumes/Maxtor/Back_up_pasteur/Kp_caen/kp_avant_2018.txt'
#working list which contains all the name of the working files
fasta_dir = '/Volumes/Maxtor/Back_up_pasteur/Kp_Caen/fasta/awked_fasta'
#directory which contains the fasta
outputdir = '/Volumes/Maxtor/Back_up_pasteur/Kp_Caen/Analyses/card_db' #path to output files
fasta_extension = 'awked.fasta'

#listing all the files which will be working on
travail = []
with open(liste, 'r') as filin:
    for nom in filin:
        travail.append(nom[:-1])

print("Voici la liste des fichiers sur lesquels vous allez travailler", travail)
print("Le nombre de fichier est de : {}".format(len(travail)))

for nom in travail:
    fasta = "{}/{}.{}".format(fasta_dir, nom, fasta_extension)
    subprocess.run(["rgi", "main", "--input_sequence", "{}".format(fasta), "--output_file",\
     "{}/{}".format(outputdir,nom), "--input_type", "contig"])


for fichier in glob.glob("{}/*.temp*".format(outputdir)):
    subprocess.run(["rm", fichier])
