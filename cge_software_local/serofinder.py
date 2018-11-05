#!/usr/bin/env python3
## Octobre 2018
### Gravey François

### module Loading
import subprocess
import re
import json
from pprint import pprint
import os.path
from shutil import copyfile
from argparse import ArgumentParser

def best_identity(dico):
    """In a dictionary, looked for the key for which the identity rate is the highest
    Input = dictionary created by parsing json file form cge script
    Output = key of the dictionary """

    id = 0
    best_key = ''
    for clef in dico.keys():
        if dico[clef]['identity'] > id:
            best_key = clef
            id = dico[clef]['identity']

    return best_key

def serofinder(liste,fasta_dir,outputdir,fasta_extension,identity,coverage):
    """ Use the serotypefinder.py maded by the cge team in order to predict the serotype of the strains
    Inputs : liste = .txt file which contains the name of the strains you wanted to work on
             fasta_dir = directory which contains the .fasta files of the strains
             outputdir = directory where you want to put the output files in
             fasta_extension = extension of the fasta file .fasta, .fsa, .agp.fasta, etc
             identity = identity rate you wanted to use
             covrage = coverage rate you wanted to use

    Outputs : """

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

    #Container definitions which will contain all the information extract from serotypefinder script
    sero = ["Souche;O;H;O identity;H identity;O/H Coverage;O Position;H Position\n"]

    # Information
    print("#########################################################")
    print("############ Launching serotypefinder ###################")
    print("#########################################################")

    #for each strain in the list, serotypefinder is executed
    for nom in travail:
        #Variable definition
        fasta = "{}{}{}".format(fasta_dir, nom, fasta_extension)

        #Information
        print("--> {}{}".format(nom, fasta_extension))

        #serotypefinder execution
        subprocess.run(["python", "{}".format(serotypefinder), "-i",\
         "{}".format(fasta), "-o", "{}".format(outputdir), "-mp", "blastn",\
         "-p", "{}".format(database), "-tmp", "{}".format(outputdir), "-l", "{}".format(coverage),\
         "-t", "{}".format(identity), "-x", "-q"])

        # Made a copy of the results.txt file form serotypefinder script
        old_file = "{}results.txt".format(outputdir)
        new_file = "{}{}_serotypefinder.txt".format(outputdir,nom)
        if os.path.isfile(old_file) == True:
            copyfile(old_file,new_file)


        #Parsing the serotypefinder results in the .json file and add all the information in the data variable
        with open('{}data.json'.format(outputdir), 'r') as filin:
            data = json.load(filin)

        ## Variables definitions
        res_O = data['serotypefinder']['results']['O_type'] #all the results obtained for O antigen
        res_H = data['serotypefinder']['results']['H_type'] #all the results obtained for H antigen

        #Working on results and checking if there is not an absence of results for each antigen
        if res_O == 'No hit found' or res_O == {}:
        #O results:
            coverage_O = "-"
            O_gene = "ND"
            O_id = "-"
            O_pos = "-"

        else:
        #O results:
            best_O = res_O[best_identity(res_O)] # looking for the best resuts thanks to best_identity funciton
            #Variables definition
            coverage_O = "{:.1f}".format(int(best_O['HSP_length'])/int(best_O['template_length'])*100)
            O_gene = "{}_{}".format(best_O['serotype'], best_O['gene'])
            O_id = best_O['identity']
            O_pos = best_O['contig_name']

        if res_H == 'No hit found' or res_H == {}:
        #H results:
            coverage_H = "-"
            H_gene = "ND"
            H_id = "-"
            H_pos = "-"

        else:
        #H results:
            best_H = res_H[best_identity(res_H)] # looking for the best resuts thanks to best_identity funciton
            #Variables definition
            coverage_H = "{:.1f}".format(int(best_H['HSP_length'])/int(best_H['template_length'])*100)
            H_gene = "{}_{}".format(best_H['serotype'], best_H['gene'])
            H_id = best_H['identity']
            H_pos = best_H['contig_name']


        #O and H antigens information appended to the sero container
        sero.append("{};{};{};{};{};{}/{};{};{}\n".format(nom, O_gene, H_gene, O_id,\
        H_id, coverage_O, coverage_H, O_pos, H_pos))

    return sero

if __name__ == "__main__":

    # Parse command line options
    ##########################################################################
    parser = ArgumentParser()
    parser.add_argument("-l", "--list", dest="list", \
    help="list which contains all the name of the strains", default='')
    parser.add_argument("-o", "--outputPath", dest="out_path",\
    help="Path to blast output", default='')
    parser.add_argument("-f", "--fasta_path", dest="fasta_path_dir",\
    help="Path to the directory which contains all the fasta to analyse", default='')
    parser.add_argument("-e", "--extension", dest="extension",\
    help="fasta extension such as .fasta .fa .awked.fasta .agp.fasta", default='')
    parser.add_argument("-filename", "--filename", dest="filename",\
    help="summary output file name", default='serotypefinder')
    parser.add_argument("-id", "--identity", dest="identity",\
    help="minimum identity rate float", default='0.8')
    parser.add_argument("-c", "--coverage", dest="coverage",\
    help="minimum coverage rate float", default='0.8')
    args = parser.parse_args()

    ###############################################################################
    # Variables difinition
    liste = args.list
    fasta_dir = args.fasta_path_dir
    outputdir = args.out_path
    fasta_extension = args.extension
    nom_fichier = args.filename
    identity = args.identity
    coverage = args.coverage
    filename = args.filename
#
    ####changing path
    # liste = '/Volumes/Maxtor/Back_up_pasteur/blse_ecoli_2018_ajout/blse_ajout.txt'
    #working list which contains all the name of the working files
    # fasta_dir = '/Volumes/Maxtor/Back_up_pasteur/blse_ecoli_2018_ajout/fasta/agp_fasta/'
    #directory which contains the fasta
    # outputdir = '/Users/Francois/Desktop/essai_sero/' #path to output files
    # fasta_extension = '.agp.fasta'
    # identity = '0.7'
    # coverage = '0.7'
    # filename = "serofinder"
#
    with open("{}{}.csv".format(outputdir,filename), "w") as filout:
        for ligne in serofinder(liste,fasta_dir,outputdir,fasta_extension,\
        identity,coverage):
            filout.write(ligne)

    print("##################################################")
    print("Cleaning process")
    subprocess.run(["rm", "-rf", "{}tmp".format(outputdir)])
    subprocess.run(["rm", "{}Serotype_allele_seq.fsa".format(outputdir)])
    subprocess.run(["rm", "{}Hit_in_genome_seq.fsa".format(outputdir)])
    subprocess.run(["rm", "{}results_tab.tsv".format(outputdir)])
    subprocess.run(["rm", "{}results.txt".format(outputdir)])
    subprocess.run(["rm", "{}data.json".format(outputdir)])
    print("####################################################")
    print("End of the process")
    print("####################################################")
