#!/usr/bin/env python3
## Septembre 2018
### Gravey François

### module Loading
import subprocess
import re

##### changing path
liste = '/Volumes/Maxtor/Back_up_pasteur/kp_caen_2018/kp_caen_2018.txt'
#working list which contains all the name of the working files
fasta_dir = '/Volumes/Maxtor/Back_up_pasteur/kp_caen_2018/fasta/agp_fasta'
#directory which contains the fasta
outputdir = '/Volumes/Maxtor/Back_up_pasteur/kp_caen_2018/analyses/pmlst' #path to output files
fasta_extension = 'agp.fasta'
#### unchanging path
pmlstfinder = '/Users/Francois/cge_softwares/pmlst/pmlst.pl'
database = '/Users/Francois/cge_data_bases/pmlst_db'
blast_path = '/Users/Francois/cge_softwares/blast-2.2.26'


#listing all the files which will be working on
travail = []
with open(liste, 'r') as filin:
    for nom in filin:
        travail.append(nom[:-1])

print("Voici la liste des fichiers sur lesquels vous allez travailler", travail)
print("Le nombre de fichier est de : {}".format(len(travail)))

for nom in travail:
    fasta = "{}/{}.{}".format(fasta_dir, nom, fasta_extension)
    outputdir_spe = "{}/{}".format(outputdir,nom)
    subprocess.run(["mkdir", outputdir_spe])
    print("##################################################")
    print("working on {}.{} file".format(nom, fasta_extension))
    print("##################################################")
    subprocess.run(["perl", "{}".format(pmlstfinder), "-d", "{}".format(database),\
    "-i", "{}".format(fasta), "-o", "{}".format(outputdir_spe), "-s", "incf",\
    "-b", "{}".format(blast_path)])
