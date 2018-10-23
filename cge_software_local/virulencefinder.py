#!/usr/bin/env python3
## Octobre 2018
### Gravey François

### module Loading
import subprocess
import re

##### changing path
liste = '/Volumes/Maxtor/lugdunensis/liste_souches_lugdu_test.txt'
#working list which contains all the name of the working files
fasta_dir = '/Volumes/Maxtor/lugdunensis/fasta_assembled/'
#directory which contains the fasta
outputdir = '/Volumes/Maxtor/lugdunensis/analyses/virulencefinder' #path to output files
fasta_extension = '.fasta'
schema = 's.aureus_exoenzyme'

#### unchanging path
virulencefinder = '/Users/Francois/cge_softwares/virulencefinder/virulencefinder.pl'
database = '/Users/Francois/cge_data_bases/virulencefinder_db'
blast_path = "/Users/Francois/cge_softwares/blast-2.2.26"

#listing all the files which will be working on
travail = []
with open(liste, 'r') as filin:
    for nom in filin:
        travail.append(nom[:-1])

print("Voici la liste des fichiers sur lesquels vous allez travailler", travail)
print("Le nombre de fichier est de : {}".format(len(travail)))

#On parcours la liste et pour chaque nom, on execute fimtyper
for nom in travail:
    fasta = "{}{}{}".format(fasta_dir, nom, fasta_extension)

    print("##################################################")
    print("working on {}{} file".format(nom, fasta_extension))
    print("##################################################")

    subprocess.run(["perl", "{}".format(virulencefinder), "-d", "{}".format(database), \
    "-i", "{}".format(fasta), "-o", "{}".format(outputdir), "-k", "85.00", \
    "-b", "{}".format(blast_path), "-s", "{}".format(schema)])
