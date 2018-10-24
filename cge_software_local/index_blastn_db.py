#!/usr/bin/env python3
## Aout 2018
### Gravey François

### module Loading
import subprocess
from glob import glob
from argparse import ArgumentParser
import os.path


def mlst_db_blastn_index(db_path):

    ## listing all the directories:
    for dir in glob("{}*/".format(db_path)):
        #Variables definitions
        nom = dir.split("/")[-2]
        fichier = "{}{}.fsa".format(dir,nom)
        # If the file exist, then using makeblastdb command for indexing
        if os.path.isfile(fichier) == True:
            subprocess.run(["makeblastdb", "-in", "{}".format(fichier), "-dbtype",\
            "nucl"])

if __name__ == "__main__":

    # PARSE COMMAND LINE OPTIONS
    ##########################################################################
    parser = ArgumentParser()
    parser.add_argument("-db", "--database", dest="database",\
    help="Indicate the path of the database", default='')
    args = parser.parse_args()

    ###############################################################################
    # Variables difinition
    db_path = args.database

    ###############################################################################
    # Launching function
    mlst_db_blastn_index(db_path)
