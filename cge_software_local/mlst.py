#!/usr/bin/env python3
## Novembre 2018
### Gravey François

### module Loading
import subprocess
import re
from argparse import ArgumentParser
import os.path
from shutil import copyfile
import json


def mlst_strain(nom,fasta_dir,fasta_extension,outputdir,specie):

    """ Use the mlst.py created by cge in order to determine the Sequence Type of the strains
    Inputs : liste = .txt file which contains the name of the strains you wanted to work on
             fasta_dir = directory which contains the .fasta files of the strains
             outputdir = directory where you want to put the output files in
             fasta_extension = extension of the fasta file .fasta, .fsa, .agp.fasta, etc
             specie = name of the specie you are working on, see the text file in the script directory

    Outputs : return a list which contains the strain name, the ST and the alleles names and numbers"""

    #unchanging paths
    mlst = '/Users/Francois/cge_softwares/mlst/mlst.py'
    mlst_db = '/Users/Francois/cge_data_bases/mlst_db'

    #Variable definition
    fasta = "{}{}{}".format(fasta_dir,nom,fasta_extension)

    #Information
    print("---> {}{}".format(nom, fasta_extension))

    #Launching mlst.py script from cge internet site
    subprocess.run(["python", "{}".format(mlst), "-i", "{}".format(fasta),\
    "-o", "{}".format(outputdir), "-s", "{}".format(specie), "-p",\
    "{}".format(mlst_db), "-t", "{}".format(outputdir), "-mp", "blastn",\
    "-x", "-q"])

    #Made a copy of the results.txt file form mlst.py script
    old_file = "{}results.txt".format(outputdir)
    new_file = "{}{}_mlst.txt".format(outputdir,nom)
    if os.path.isfile(old_file) == True:
        copyfile(old_file,new_file)


    #Parsing the serotypefinder results in the .json file and add all the information in the data variable
    with open('{}data.json'.format(outputdir), 'r') as filin:
        data = json.load(filin)

    #Variables definition
    plus_proche_st = []
    ST = data["mlst"]["results"]["sequence_type"]

    # if ST is Unknown adding the nearest ST when is exists
    if ST == "Unknown":
        plus_proche_st.append(data["mlst"]["results"]["nearest_sts"])
    else:
        plus_proche_st.append("-")

    if not plus_proche_st:
        plus_proche_st.append("-")

    #Looking for allele results
    allele = []
    for clef in data["mlst"]["results"]["allele_profile"]:
        allele.append(data["mlst"]["results"]["allele_profile"][clef]["allele_name"].replace("_", "-"))

    #Creating the res container for each strain
    res = "{};ST{};{};{}\n".format(nom,ST,",".join(sorted(allele)),"".join(plus_proche_st))

    #End of the function
    return(res)

def mlst_all(liste,fasta_dir,fasta_extension,outputdir,specie):

    """ Use the mlst.py created by cge in order to determine the Sequence Type of the strains
    Inputs : liste = .txt file which contains the name of the strains you wanted to work on
             fasta_dir = directory which contains the .fasta files of the strains
             outputdir = directory where you want to put the output files in
             fasta_extension = extension of the fasta file .fasta, .fsa, .agp.fasta, etc
             specie = name of the specie you are working on, see the text file in the script directory

    Outputs : for each strain, the results of mlst.py are collected into a
              dedicated file which is named as strain_name_mlst.txt
              a summary list which contains all the results for all the tested strains"""

    #listing all the files which will be working on
    travail = []
    with open(liste, 'r') as filin:
        for nom in filin:
            travail.append(nom[:-1])

    print("Voici la liste des fichiers sur lesquels vous allez travailler", travail)
    print("Le nombre de fichier est de : {}".format(len(travail)))

    #Container definitions which will contain all the information extract from serotypefinder script
    mlst = ["Souche;ST;alleles;plus_proche_st\n"]

    # Information
    print("#########################################################")
    print("############ Launching mlstfinder #######################")
    print("#########################################################")

    #Launching mlst_strain function
    for nom in travail:
        mlst.append(mlst_strain(nom,fasta_dir,fasta_extension,outputdir,specie))
#
    #Cleaning process
    # subprocess.run(["rm", "-rf", "{}tmp".format(outputdir)])
    # subprocess.run(["rm", "{}MLST_allele_seq.fsa".format(outputdir)])
    # subprocess.run(["rm", "{}Hit_in_genome_seq.fsa".format(outputdir)])
    # subprocess.run(["rm", "{}results_tab.tsv".format(outputdir)])
    # subprocess.run(["rm", "{}results.txt".format(outputdir)])
    # subprocess.run(["rm", "{}data.json".format(outputdir)])

    #End of the function
    return(mlst)

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
    help="Any informations such as .agp.fasta .awked.fa .scfd.fasta .fasta .fa",\
     default='')
    parser.add_argument("-o", "--outputPath", dest="out_path",help="Path to blast output",\
     default='')
    parser.add_argument("-s", "--specie", dest="specie",help="Which specie are you working on ?",\
     default='')
    parser.add_argument("-filename", "--filename", dest="filename",help="name of the summary file",\
     default='mlst_finder')
    args = parser.parse_args()

    #Defining varibales

    liste = args.list
    fasta_dir = args.fasta_dir
    outputdir = args.out_path
    fasta_extension = args.extension
    specie = args.specie
    filename = args.filename

    # Launching mlst_all function
    with open("{}{}.csv".format(outputdir,filename), "w") as filout:
        for ligne in mlst_all(liste,fasta_dir,fasta_extension,outputdir,specie):
            filout.write(ligne)
