#!/usr/bin/env python3
## Novembre 2018
### Gravey François

import subprocess

def specie(nom,fasta_dir,fasta_extension):
    #Variable definition
    fasta = "{}{}{}".format(fasta_dir,nom,fasta_extension)

    #Unvariable path
    rmlst = '/Users/Francois/Programs/PubMLST/rMLST.py'

    #Lauching the rmlst script from PubMLST
    res = subprocess.check_output(["python", "{}".format(rmlst), "--file", "{}".format(fasta)])
    res = res.decode("utf-8")
    return("{};{};{};{};{}\n".format(nom, res.split('\n')[0].split(' ')[1],\
     res.split('\n')[1].split(":")[1],res.split('\n')[2].split(":")[1],\
     res.split('\n')[3]))

def specie_all(liste,fasta_dir,fasta_extension):
    #Container creation : list which will contain all the strains name of the project
    travail = []
    #Reading the name of the strains listed into a .txt file
    with open(liste, "r") as filin:
        for nom in filin:
            travail.append(nom[:-1])

    #Information
    print("\n")
    print("Voici les souches sur lesquelles vous allez travailler {}\n".format(travail))
    print("Le nombre de souches est de {}\n".format(len(travail)))
    print("####################################################################")
    print("############## launch the rMLST script #############################")
    print("####################################################################")

    # for each strain lauching the specie function
    with open("{}{}.csv".format(outputdir,filename), "w") as filout:
        filout.write("Souche;Rank;Taxon;Support;Taxonomy\n")
        for nom in travail:
            # Information
            print("---> {}{}".format(nom, fasta_extension))

            #launch the function
            filout.write(specie(nom,fasta_dir,fasta_extension))

if __name__ == "__main__":

    # Parse command line options
    ##########################################################################
    parser = ArgumentParser()
    parser.add_argument("-l", "--list", dest="list", \
    help="list which contains all the name of the strains", default='')
    parser.add_argument("-f", "--fasta_path", dest="fasta_path_dir",\
    help="Path to the directory which contains all the fasta to analyse", default='')
    parser.add_argument("-e", "--extension", dest="extension",\
    help="fasta extension such as .fasta .fa .awked.fasta .agp.fasta", default='')
    parser.add_argument("-o", "--outputPath", dest="out_path",\
    help="Path to blast output", default='')
    parser.add_argument("-filename", "--filename", dest="filename",\
    help="summary output file name", default='blastn_best')
    args = parser.parse_args()

    ###############################################################################
    # Variables difinition
    liste = args.list
    fasta_dir = args.fasta_path_dir
    fasta_extension = args.extension
    outputdir = args.out_path
    filename = args.filename

nom = 'P1-01'
liste = "/Users/Francois/Documents/projets/ecoli/Rea_neonat/rea_test.txt"
fasta_dir = '/Users/Francois/Documents/projets/ecoli/Rea_neonat/fasta/awked_fasta/'
fasta_extension = '.awked.fasta'
outputdir = '/Users/Francois/Desktop/essai/'
filename = 'essai_specie'
specie_all(liste,fasta_dir,fasta_extension)
