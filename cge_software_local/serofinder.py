#!/usr/bin/env python3
## Aout 2018
### Gravey François

### module Loading
import subprocess
import re
import json
from pprint import pprint

def best_identity(dico):
    id = 0
    best_key = ''
    for clef in dico.keys():
        if dico[clef]['identity'] > id:
            best_key = clef
            id = dico[clef]['identity']

    return best_key

def serofinder(liste,fasta_dir,outputdir,fasta_extension,identity,coverage):
    #### unchanging path
    serotypefinder = '/Users/Francois/cge_softwares/serotypefinder/serotypefinder.py'
    database = '/Users/Francois/cge_data_bases/serotypefinder_db'


    #listing all the files which will be working on
    travail = []
    with open(liste, 'r') as filin:
        for nom in filin:
            travail.append(nom[:-1])

    print("Voici la liste des fichiers sur lesquels vous allez travailler", travail)
    print("Le nombre de fichier est de : {}".format(len(travail)))

    sero = ["Souche;O;H;O identity;H identity;O/H Coverage;O Position;H Position\n"]
    #On parcours la liste et pour chaque nom, on execute fimtyper
    for nom in travail:
        fasta = "{}{}{}".format(fasta_dir, nom, fasta_extension)
        print("##################################################")
        print("working on {}{} file".format(nom, fasta_extension))
        print("##################################################")
        subprocess.run(["python", "{}".format(serotypefinder), "-i",\
         "{}".format(fasta), "-o", "{}".format(outputdir), "-mp", "blastn",\
         "-p", "{}".format(database), "-tmp", "{}".format(outputdir), "-l", "{}".format(coverage),\
         "-t", "{}".format(identity), "-x", "-q"])

        with open('{}data.json'.format(outputdir), 'r') as filin:
            data = json.load(filin)

        ## Variables definitions
        res_O = data['serotypefinder']['results']['O_type']
        res_H = data['serotypefinder']['results']['H_type']

        if res_O == 'No hit found':
        #O results:
            coverage_O = "-"
            O_gene = "ND"
            O_id = "-"
            O_pos = "-"

        else:
        #O results:
            best_O = res_O[best_identity(res_O)]
            coverage_O = "{:.1f}".format(int(best_O['HSP_length'])/int(best_O['template_length'])*100)
            O_gene = "{}_{}".format(best_O['serotype'], best_O['gene'])
            O_id = best_O['identity']
            O_pos = best_O['contig_name']

        if res_H == 'No hit found':
        #H results:
            coverage_H = "-"
            H_gene = "ND"
            H_id = "-"
            H_pos = "-"

        else:
        #H results:
            best_H = res_H[best_identity(res_H)]
            coverage_H = "{:.1f}".format(int(best_H['HSP_length'])/int(best_H['template_length'])*100)
            H_gene = "{}_{}".format(best_H['serotype'], best_H['gene'])
            H_id = best_H['identity']
            H_pos = best_H['contig_name']


        #Information append to sero
        sero.append("{};{};{};{};{};{}/{};{};{}\n".format(nom, O_gene, H_gene, O_id,\
        H_id, coverage_O, coverage_H, O_pos, H_pos))


##### changing path
liste = '/Volumes/Maxtor/Back_up_pasteur/blse_ecoli_2018_ajout/blse_ajout.txt'
#working list which contains all the name of the working files
fasta_dir = '/Volumes/Maxtor/Back_up_pasteur/blse_ecoli_2018_ajout/fasta/agp_fasta/'
#directory which contains the fasta
outputdir = '/Users/Francois/Desktop/essai_sero/' #path to output files
fasta_extension = '.agp.fasta'
identity = '0.9'
coverage = '0.9'
serofinder(liste,fasta_dir,outputdir,fasta_extension,identity,coverage)
#
###unchanging path
# serotypefinder = '/Users/Francois/cge_softwares/serotypefinder/serotypefinder.py'
# database = '/Users/Francois/cge_data_bases/serotypefinder_db'
#
#
#listing all the files which will be working on
# travail = []
# with open(liste, 'r') as filin:
    # for nom in filin:
        # travail.append(nom[:-1])
#
# print("Voici la liste des fichiers sur lesquels vous allez travailler", travail)
# print("Le nombre de fichier est de : {}".format(len(travail)))
#
# sero = ["Souche;O;H;O identity;H identity;O/H Coverage;O Position;H Position\n"]
#On parcours la liste et pour chaque nom, on execute fimtyper
# for nom in travail:
    # fasta = "{}{}{}".format(fasta_dir, nom, fasta_extension)
    # print("##################################################")
    # print("working on {}{} file".format(nom, fasta_extension))
    # print("##################################################")
    # subprocess.run(["python", "{}".format(serotypefinder), "-i",\
     # "{}".format(fasta), "-o", "{}".format(outputdir), "-mp", "blastn",\
     # "-p", "{}".format(database), "-tmp", "{}".format(outputdir), "-l", "0.9",\
     # "-t", "0.9", "-x"])
#
    # with open('{}data.json'.format(outputdir), 'r') as filin:
        # data = json.load(filin)
#
    #Variables definitions
    # res_O = data['serotypefinder']['results']['O_type']
    # res_H = data['serotypefinder']['results']['H_type']
#
    # if res_O == 'No hit found':
    #O results:
        # coverage_O = "-"
        # O_gene = "ND"
        # O_id = "-"
        # O_pos = "-"
#
    # else:
    #O results:
        # best_O = res_O[best_identity(res_O)]
        # coverage_O = "{:.1f}".format(int(best_O['HSP_length'])/int(best_O['template_length'])*100)
        # O_gene = "{}_{}".format(best_O['serotype'], best_O['gene'])
        # O_id = best_O['identity']
        # O_pos = best_O['contig_name']
#
    # if res_H == 'No hit found':
    #H results:
        # coverage_H = "-"
        # H_gene = "ND"
        # H_id = "-"
        # H_pos = "-"
#
    # else:
    #H results:
        # best_H = res_H[best_identity(res_H)]
        # coverage_H = "{:.1f}".format(int(best_H['HSP_length'])/int(best_H['template_length'])*100)
        # H_gene = "{}_{}".format(best_H['serotype'], best_H['gene'])
        # H_id = best_H['identity']
        # H_pos = best_H['contig_name']
#
#
    #Information append to sero
    # sero.append("{};{};{};{};{};{}/{};{};{}\n".format(nom, O_gene, H_gene, O_id,\
    # H_id, coverage_O, coverage_H, O_pos, H_pos))
#
#
    # print("##################################################")
    # print("Cleaning process")
    # subprocess.run(["rm", "-rf", "{}/tmp".format(outputdir)])
    # subprocess.run(["rm", "{}/Serotype_gene_seq.fsa".format(outputdir)])
    # subprocess.run(["rm", "{}/Hit_in_genome_seq.fsa".format(outputdir)])
    # subprocess.run(["rm", "{}/results_tab.txt".format(outputdir)])
    # subprocess.run(["rm", "{}/results.txt".format(outputdir)])
    # subprocess.run(["rm", "{}/results_table.txt".format(outputdir)])
    # print("End of the process")

#
# with open("{}/serotype_all_souches.csv".format(outputdir), 'w') as filout:
    # for l in sero:
        # filout.write(l)
        #
