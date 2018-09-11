#!/usr/bin/env python3
## Aout 2018
### Gravey François

### module Loading
import subprocess
import re

##### changing path
liste = '/Volumes/shigella-ngs/EcCaen/Rea_neonat/fichiers.txt'
#working list which contains all the name of the working files
fasta_dir = '/Volumes/shigella-ngs/EcCaen/Rea_neonat/fasta/awked_fasta'
#directory which contains the fasta
outputdir = '/Users/Francois/Desktop/test_serotype' #path to output files
fasta_extension = 'awked.fasta'

#### unchanging path
serotyper = '/Users/Francois/cge_softwares/serotypefinder/serotypefinder.pl'
database = '/Users/Francois/cge_data_bases/serotypefinder_db'
blast_path = "/Users/Francois/cge_softwares/blast-2.2.26"

#listing all the files which will be working on
travail = []
with open(liste, 'r') as filin:
    for nom in filin:
        travail.append(nom[:-1])

print("Voici la liste des fichiers sur lesquels vous allez travailler", travail)
print("Le nombre de fichier est de : {}".format(len(travail)))
header = "souche\tSerotype predit\n"
sero = []
#On parcours la liste et pour chaque nom, on execute fimtyper
for nom in travail:
    fasta = "{}/{}.{}".format(fasta_dir, nom, fasta_extension)
    print("##################################################")
    print("working on {}.{} file".format(nom, fasta_extension))
    print("##################################################")
    subprocess.run(["perl", "{}".format(serotyper), "-d", "{}".format(database), \
    "-i", "{}".format(fasta), "-o", "{}".format(outputdir), "-k", "85.00", "-l",\
     "0.60", "-b", "{}".format(blast_path), "-s", "ecoli"])

    resultat = []
    with open('{}/results_table.txt'.format(outputdir), 'r') as filin:
        with open("{}/serotype_{}.txt".format(outputdir,nom), 'w') as filout:
            for ligne in filin:
                filout.write(ligne)
                resultat.append(ligne[:-1])

    regex = re.compile("No serotype predicted\.")

    if regex.search(resultat[5]):
        #print("{}_({}%:{})".format(resultat[2].split('\t')[5],\
        #resultat[2].split('\t')[1],resultat[2].split('\t')[2]))
        sero.append("{}\tO_unpredicted:{}_({}%:{})\n".format(nom,resultat[2].split('\t')[5],\
        resultat[2].split('\t')[1],resultat[2].split('\t')[2]))
    else:
        if resultat[7]:
            sero.append("{}\t{}_({} {}%:{}, {} {}%:{}):{}_({}%:{})\n".format(nom,\
            resultat[6].split('\t')[5],resultat[6].split('\t')[0],\
            resultat[6].split('\t')[1],resultat[6].split('\t')[2],\
            resultat[7].split('\t')[0],resultat[7].split('\t')[1],\
            resultat[7].split('\t')[2],resultat[2].split('\t')[5],\
            resultat[2].split('\t')[1],resultat[2].split('\t')[2]))
        else:
            sero.append("{}\t{}_({} {}%:{}):{}_({}%:{})\n".format(nom,\
            resultat[6].split('\t')[5],resultat[6].split('\t')[0],\
            resultat[6].split('\t')[1],resultat[6].split('\t')[2],\
            resultat[2].split('\t')[5],resultat[2].split('\t')[1],\
            resultat[2].split('\t')[2]))

    print("##################################################")
    print("Cleaning process")
    subprocess.run(["rm", "-rf", "{}/tmp".format(outputdir)])
    subprocess.run(["rm", "{}/Serotype_gene_seq.fsa".format(outputdir)])
    subprocess.run(["rm", "{}/Hit_in_genome_seq.fsa".format(outputdir)])
    subprocess.run(["rm", "{}/results_tab.txt".format(outputdir)])
    subprocess.run(["rm", "{}/results.txt".format(outputdir)])
    subprocess.run(["rm", "{}/results_table.txt".format(outputdir)])
    print("End of the process")
    print("##################################################")

with open("{}/serotype_all_souches.txt".format(outputdir), 'w') as filout:
    filout.write(header)
    for l in sero:
        filout.write(l)
