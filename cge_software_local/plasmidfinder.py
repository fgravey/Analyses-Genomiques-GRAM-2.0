#!/usr/bin/env python3
## Septembre 2018
### Gravey François

### module Loading
import subprocess
import re
from argparse import ArgumentParser

def plasmidfinder_all(liste,fasta_dir, fasta_extension, outputdir, specie):

    """ Looking for plasmids carried by strains using the plasmidfinder.pl script
    - Inputs : liste = .txt file which contains all the studied strains
               fasta_dir = path to the directory which contains all the fasta files
               fasta_extension = extension of the fasta files such as .fasta or .fna
               outputdir = path to te directory which will contain all the results
    -Outputs : for each strain, the results of plasmidfinder.pl are collected into a
               dedicated directory which is named as the strain name
               a summary file which contains all the results for all the tested strains
               """

    #### unchanging path
    plasmidfinder = '/Users/Francois/cge_softwares/plasmidfinder/plasmidfinder.py'
    database = '/Users/Francois/cge_data_bases/plasmidfinder_db'

    #listing all the files which will be working on
    travail = []
    with open(liste, 'r') as filin:
        for nom in filin:
            travail.append(nom[:-1])

    # Informations
    print("Voici la liste des fichiers sur lesquels vous allez travailler", travail)
    print("Le nombre de fichier est de : {}".format(len(travail)))

    print("##################################################")
    print("####### Launching plasmidfinder Script ###########")
    print("##################################################")

    # Creation of a container which will contain all the results
    resultats = ["Strain;Database;Plasmide;Identities;Alignment.Length;Template.Length;\
    Position.in.reference;Contig;Position.in.contig;Note;Accession no.\n"]

    for nom in travail:
        #Variables creation
        fasta = "{}{}{}".format(fasta_dir, nom, fasta_extension)
        outputdir2 = "{}{}".format(outputdir,nom)

        #Make a specific directory form each strain tested
        subprocess.run(["mkdir", outputdir2])

        #Information
        print("---> {}{}".format(nom, fasta_extension))

        #Launching the plasmidfinder.pl script
        subprocess.run(["python", "{}".format(plasmidfinder), "-p", "{}".format(database),\
         "-i", "{}".format(fasta), "-o", "{}".format(outputdir2, "-d", "{}".format(specie))])

        #Creation of a container which will contains the results for each strain
        plasmide = []

        #Parsing the output of the plasmidfinder.pl script and extract somes relevant results
        with open("{}/results_tab.txt".format(outputdir2),"r") as filin:
            for ligne in filin:
                plasmide.append(ligne[:-1])
            if len(plasmide) == 1:
                resultats.append("{};{};no plasmid found;-;-;-;-;-;-;-;-\n"\
                .format(nom,"Enterobacteriaceae"))
            else:
                for i in range(1,len(plasmide)):
                    resultats.append("{};{}\n".format(nom,plasmide[i].replace("\t",";")))

    #End of the function
    return(resultats)

if __name__ == "__main__":
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
    parser.add_argument("-s", "--specie", dest="specie",\
    help="Which type of database would you lie to use ? enterobacteriaceae or gram_positive", default='')
    parser.add_argument("-o", "--outputPath", dest="out_path",help="Path to blast output", default='')
    args = parser.parse_args()

    # Defining varibales

    liste = args.list
    fasta_dir = args.fasta_dir
    outputdir = args.out_path
    fasta_extension = args.extension
    specie = args.specie

    # Using plasmidfinder_all function
    with open("{}/plasmidfinder_results.csv".format(outputdir), "w") as filout:
        for res in plasmidfinder_all(liste,fasta_dir, fasta_extension, outputdir,specie):
            filout.write(res)
