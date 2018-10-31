#!/usr/bin/env python3
## Octobre 2018
### Gravey François

### module Loading
import subprocess
import re
from argparse import ArgumentParser
import glob
import os

def phylogroup(nom,fasta_dir,fasta_extension,database,outputdir):
    #Varibales definitions
    fasta = "{}{}{}".format(fasta_dir,nom,fasta_extension)
    tmp = "{}tmp".format(outputdir)

    #Constant paths
    query = "/Users/Francois/Programs/ClermonTyping/data/primers.fasta"
    clermont = "/Users/Francois/Programs/ClermonTyping/bin/clermont.py"

    #Execution of several steps :
    ## First creation of the database form the fasta file in the temp directory
    subprocess.run(["makeblastdb", "-in", "{}".format(fasta),"-input_type", \
    "fasta", "-out", "{}/{}".format(tmp,nom) ,"-dbtype", "nucl"])

    ## Execution of blastn between the created database and the clermont query
    subprocess.run(["blastn", "-query", "{}".format(query), "-perc_identity", \
    "90" ,"-task", "blastn", "-outfmt", "5" ,"-db", "{}/{}".format(tmp,nom), \
    "-out", "{}/{}.xml".format(outputdir,nom)])

    ## Interpretation of the results using the clermont.py script
    res = subprocess.check_output(["/Users/Francois/Programs/ClermonTyping/bin/clermont.py",\
     "-x", "{}/{}.xml".format(outputdir,nom)])
    ### res is a variable which contains the stdout of the clermont.py script
    res = res.decode("utf-8")
    res = "{};{}".format(nom, res.replace("\t", ";"))

    #End of the function
    return(res)

def phylogroup_list(liste,fasta_dir,fasta_extension,database,outputdir):
    # listing all the strains of the project
    travail = []
    with open(liste, 'r') as filin:
        for nom in filin:
            travail.append(nom[:-1])

    #Information
    print("Voici les souches sur lesquelles vous allez travailler {}".format(travail))
    print("Leur nombre est de {}".format(len(travail)))
    print("########################################################")
    print("################# Phylogroup Execution #################")
    print("########################################################")
    print("\n")

    #For each strain, using the phylogroup function and add the results into resultats
    resultats = ["Souche;Genes;Resultats genes;gene de typage;Phylogroup\n"]
    for nom in travail:
        #Information
        print("---> {}\n".format(nom))
        resultats.append(phylogroup(nom,fasta_dir,fasta_extension,database,outputdir))

    #Cleaning process:
    subprocess.run(["rm", "-rf", "{}tmp".format(outputdir)])

    #End of the function
    return(resultats)


liste = '/Volumes/Maxtor/Back_up_pasteur/Rea_neonat/rea.txt'
outputdir = '/Users/Francois/Desktop/essai_phylogroupes/'
fasta_dir = '/Volumes/Maxtor/Back_up_pasteur/Rea_neonat/fasta/agp_fasta/'
fasta_extension = '.agp.fasta'
database = '/Users/Francois/blast_data_base/phylogroups/phylogroups.fasta'
with open("{}phylogroup.csv".format(outputdir), "w") as filout:
    for ligne in phylogroup_list(liste,fasta_dir,fasta_extension,database,outputdir):
        filout.write(ligne)
