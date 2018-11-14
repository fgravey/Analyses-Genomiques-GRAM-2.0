#!/usr/bin/env python3
#Octobre 2018
#Gravey François

import subprocess
import os
from argparse import ArgumentParser
from shutil import copyfile

##loading modules
subprocess.run(["module load SPAdes/3.12.0"], shell = True)

def spades(liste,fastq_dir,fastq_extension,outputdir):
    ## Variables difinitions:
    travail = [] #List which will contain all the name of the strains

    #Reading the list file and adding all the strain name into the travail list
    with open(liste,'r') as filin:
        for ligne in filin:
            travail.append(ligne[:-1])

    ## Information about which strains the program will working on :
    print("########################################################")
    print("############ Voici les souches sur lesquelles vous allez travailler #########################")
    print(travail)
    print("############# Leur nombre est de {}########################".format(len(travail)))
    print("\n")

    ## For each strain, launching the SPAdes program:
    for nom in travail:
        print("########################################################")
        print("################# Working on {} ##########################".format(nom))
        print("########################################################")
        print("\n")

        ### Variables definition:
        R1 = "{}{}{}".format(fastq_dir, nom, fastq_extension[0])
        R2 = "{}{}{}".format(fastq_dir, nom, fastq_extension[1])
        outputdir_assembly = "{}{}".format(outputdir,nom)
        outputdir_assembly_temp = "{}/temp".format(outputdir_assembly)

        log = "{}/log.txt".format(outputdir_assembly)
        err = "{}/err.txt".format(outputdir_assembly)


        ### Directory Creation:
        subprocess.run(["mkdir", "-p", outputdir_assembly_temp])

        #SPAdes sbatch
        subprocess.run(["sbatch" , "-J SPAdes", "--mem=65000","-o", log, "-e", err, \
        "--wrap=spades.py -1 {} -2 {} -o {} --tmp-dir {} -k {} --careful --cov-cutoff 10.0".\
        format(R1, R2, outputdir_assembly, outputdir_assembly_temp, kmer)])

def rename_scaffolds(liste,outputdir):
    ## Variables difinitions:
    travail = [] #List which will contain all the name of the strains

    #Reading the list file and adding all the strain name into the travail list
    with open(liste,'r') as filin:
        for ligne in filin:
            travail.append(ligne[:-1])

    #Rename each scaffolds.fasta like name.scfd.fasta
    for nom in travail:
        directory = "{}{}".format(outputdir,nom)
        old_file = os.path.join(directory, "scaffolds.fasta")
        new_file = os.path.join(directory, "{}.scfd.fasta".format(nom))
        os.rename(old_file, new_file)

def mv_fasta(liste,outputdir):
    ## Variables difinitions:
    travail = [] #List which will contain all the name of the strains

    #Reading the list file and adding all the strain name into the travail list
    with open(liste,'r') as filin:
        for ligne in filin:
            travail.append(ligne[:-1])

    #create directory which will contains all the assembled fasta
    subprocess.run(["mkdir", "{}fasta_assembled".format(outputdir)])

    #Rename each scaffolds.fasta like name.scfd.fasta
    for nom in travail:
        old_file = "{}{}/{}.scfd.fasta".format(outputdir,nom,nom)
        new_file = "{}fasta_assembled/{}.scfd.fasta".format(outputdir,nom,nom)
        copyfile(old_file, new_file)

if __name__ == "__main__":

    ## looking for scripts arguments
    parser = ArgumentParser()
    parser.add_argument("-l", "--list", dest="list", \
    help="list which contains all the name of the strains", default='')
    parser.add_argument("-o", "--outputdir", dest="outputdir",\
    help="Path to the directory which will contain the genomes assembly output", default='')
    parser.add_argument("-f", "--fastq_path", dest="fastq_path_dir",\
    help="Path to the directory which contains all the fastq to assemble", default='')
    parser.add_argument("-e", "--extension", dest="extension",\
    help="fastq extension such as _R1_fastq _R2_fastq", default=["_R1.fastq","_R2.fastq"])
    parser.add_argument("-filename", "--filename", dest="filename",\
    help="summary output file name", default='summary_blast')
    parser.add_argument("-k", "--kmer", dest="kmer",\
    help="list of kmer size", default='21,33,55,77')
    args = parser.parse_args()

    ## Varibales definition
    liste = args.list
    fastq_dir = args.fastq_path_dir
    outputdir = args.outputdir
    fastq_extension = args.extension
    nom_fichier = args.filename
    kmer = args.kmer


    spades(liste,fastq_dir,fastq_extension,outputdir)
    #rename_scaffolds(liste,outputdir)
    #mv_fasta(liste,outputdir)
