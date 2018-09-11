#!/usr/bin/env python3
## Aout 2018
### Gravey François

### module Loading
import subprocess

##### changing path
liste = '/Volumes/shigella-ngs-1/EcCaen/Rea_neonat/Enterobase/fasta\
/travailles_phylogenie/fasta_cluster_droit.txt'
#working list which contains all the name of the working files
fasta_dir = '/Volumes/shigella-ngs-1/EcCaen/Rea_neonat/Enterobase/fasta/travailles_phylogenie/cluster_droit'
#directory which contains the fasta
outputdir = '/Volumes/shigella-ngs-1/EcCaen/Rea_neonat/Analyses/Enterobase/Cluster_droite/fim' #path to output files
fasta_extension = 'fasta'

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
    fasta = "{}/{}.{}".format(fasta_dir, nom, fasta_extension)
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
            filout.write(resultat[3] + '\n')
            filout.write(resultat[4] + '\n')
    else:
        with open('{}/resulats_fimH_typing.txt'.format(outputdir), 'a') as filout:
            filout.write(resultat[4] + '\n')

print("##################################################")
print("Cleaning process")
subprocess.run(["rm", "-rf", "{}/tmp".format(outputdir)])
subprocess.run(["rm", "{}/Fimtype_seq.fsa".format(outputdir)])
subprocess.run(["rm", "{}/Hit_in_genome_seq.fsa".format(outputdir)])
subprocess.run(["rm", "{}/results_tab.txt".format(outputdir)])
subprocess.run(["rm", "{}/results.txt".format(outputdir)])
print("End of the process")
print("##################################################")
