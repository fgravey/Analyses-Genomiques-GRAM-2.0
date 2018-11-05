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
    ST = data["mlst"]["results"]["sequence_type"]

    #Looking for allele results
    allele = []
    for clef in data["mlst"]["results"]["allele_profile"]:
        allele.append(data["mlst"]["results"]["allele_profile"][clef]["allele_name"].replace("_", "-"))

    #Creating the res container for each strain
    res = "{};ST{};{}".format(nom,ST,",".join(allele))

    #End of the function
    return(res)


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
    parser.add_argument("-o", "--outputPath", dest="out_path",help="Path to blast output", default='')
    parser.add_argument("-s", "--specie", dest="specie",help="Which specie are you working on ?", default='')
    args = parser.parse_args()
#
    #Defining varibales
#
    # liste = args.list
    # fasta_dir = args.fasta_dir
    # outputdir = args.out_path
    # fasta_extension = args.extension
    # specie = args.specie

    nom = 'P1-01'
    fasta_dir = '/Users/Francois/Documents/projets/ecoli/Rea_neonat/fasta/awked_fasta/'
    outputdir = '/Users/Francois/Desktop/mlst/'
    fasta_extension = '.awked.fasta'
    specie = "ecoli"
    mlst_strain(nom,fasta_dir,fasta_extension,outputdir,specie)
