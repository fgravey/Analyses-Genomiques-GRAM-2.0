#!/usr/bin/env python3
## Octobre 2018
### Gravey François

#### Libraries importation
from shutil import copyfile
from argparse import ArgumentParser
from glob import glob

outputdir = "/Users/Francois/Documents/projets/ecloacae/fasta_patric/"
inputdir = "/Users/Francois/Documents/projets/ecloacae/patric_db"

for dir in glob("{}/*/".format(inputdir)):
    for file in 
