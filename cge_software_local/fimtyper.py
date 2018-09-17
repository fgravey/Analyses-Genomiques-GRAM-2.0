#!/usr/bin/env python3
## Aout 2018
### Gravey François

### module Loading
import subprocess
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

#### unchanging path
fimtyper = '/Users/Francois/cge_softwares/fimtyper/fimtyper.pl'
database = '/Users/Francois/cge_softwares/fimtyper/fimtyper_db'

#listing all the files which will be working on
travail = []
with open(liste, 'r') as filin:
    for nom in filin:
        travail.append(nom[:-1])

print("Voici la liste des fichiers sur lesquels vous allez travailler", travail)
print("Le nombre de fichier est de : {}".format(len(travail)))

#On parcours la liste et pour chaque nom, on execute fimtyper
for nom in travail:
    fasta = "{}/{}{}".format(fasta_dir, nom, fasta_extension)
    print("##################################################")
    print("working on {}.{} file".format(nom, fasta_extension))
    print("##################################################")
    subprocess.run(["perl", "{}".format(fimtyper), "-d", "{}".format(database), \
    "-i", "{}".format(fasta), "-o", "{}".format(outputdir), "-k", "95.00", "-l",\
     "0.60"])

    resultat = []
    with open('{}/results_tab.txt'.format(outputdir), 'r') as filin:
        for ligne in filin:
            resultat.append(ligne[:-1])

    if travail.index("{}".format(nom)) == 0:
        with open('{}/resulats_fimH_typing.txt'.format(outputdir), 'w') as filout:
            filout.write("{}\t{}\n".format(nom,resultat[3]))
            filout.write("{}\t{}\n".format(nom,resultat[4]))
    else:
        with open('{}/resulats_fimH_typing.txt'.format(outputdir), 'a') as filout:
            filout.write("{}\t{}\n".format(nom,resultat[4]))

print("##################################################")
print("Cleaning process")
subprocess.run(["rm", "-rf", "{}/tmp".format(outputdir)])
subprocess.run(["rm", "{}/Fimtype_seq.fsa".format(outputdir)])
subprocess.run(["rm", "{}/Hit_in_genome_seq.fsa".format(outputdir)])
subprocess.run(["rm", "{}/results_tab.txt".format(outputdir)])
subprocess.run(["rm", "{}/results.txt".format(outputdir)])
print("End of the process")
print("##################################################")
