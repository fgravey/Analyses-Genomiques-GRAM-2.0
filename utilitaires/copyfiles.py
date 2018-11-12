#!/usr/bin/env python3
## Aout 2018
### Gravey François

#### Libraries importation
from shutil import copyfile
from argparse import ArgumentParser

def copy_files(liste,inputdir,outputdir, extension):
    travail = []
    with open(liste, 'r') as filin:
        for ligne in filin:
            travail.append(ligne[:-1])

    for nom in travail:
        old_file = "{}{}{}".format(inputdir,nom,extension)
        new_file = "{}{}{}".format(outputdir,nom,extension)
        print("{} ---> {}".format(old_file, new_file))
        copyfile(old_file,new_file)

if __name__ == "__main__":

    ####################################################################################
    ################################# Main #############################################
    ####################################################################################
    # PARSE COMMAND LINE OPTIONS
    ##########################################################################
    parser = ArgumentParser()
    parser.add_argument("-l", "--list", dest="list", \
    help="list which contains all the name of the strains", default='')
    parser.add_argument("-o", "--outputPath", dest="out_path",\
    help="Path where you wanted to copy your files output", default='')
    parser.add_argument("-i", "--input_path", dest="input_path_dir",\
    help="Path to the directory which contains the files you wanted to move", default='')
    parser.add_argument("-e", "--extension", dest="extension",\
    help="extension of the files you wanted to move .fasta, .py, .fastq,...", default='')

    args = parser.parse_args()

    ###############################################################################


    # Variables difinition
    liste = args.list
    inputdir = args.input_path_dir
    outputdir = args.out_path
    extension = args.extension

    # Function execution
    copy_files(liste,inputdir,outputdir, extension)
