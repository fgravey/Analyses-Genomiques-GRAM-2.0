#!/usr/bin/env python3
## Octobre 2018
### Gravey François

### module Loading
import subprocess
import re
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from argparse import ArgumentParser
import glob
import os

nom = 'P1-01'
outputdir = '/Users/Francois/Desktop/essai_phylogroupes/'
fasta_dir = '/Volumes/Maxtor/Back_up_pasteur/Rea_neonat/fasta/agp_fasta/'
fasta_extension = '.agp.fasta'
database = '/Users/Francois/blast_data_base/phylogroups/phylogroups.fasta'
#blastn(nom,outputdir,fasta_dir,fasta_extension,database)
fasta = "{}{}{}".format(fasta_dir,nom,fasta_extension)
tmp = "{}tmp".format(outputdir)
query = "/Users/Francois/Programs/ClermonTyping/data/primers.fasta"

subprocess.run(["makeblastdb", "-in", "{}".format(fasta),"-input_type", "fasta", "-out", "{}/{}".format(tmp,nom) ,"-dbtype", "nucl"])
subprocess.run(["blastn", "-query", "{}".format(query), "-perc_identity", "90" ,"-task",\
 "blastn", "-outfmt", "5" ,"-db", "{}/{}".format(tmp,nom), "-out", "{}/{}.xml".format(outputdir,nom)])
phylogroup = ''
subprocess.run(["/Users/Francois/Programs/ClermonTyping/bin/clermont.py", "-x", "{}/{}.xml".format(outputdir,nom)],stdout = phylogroup)
# with open("{}{}_phylogroups.txt".format(outputdir,nom), "w") as filout:
    # sortie = "{}\t{}".format(nom,)
    # subprocess.run(["/Users/Francois/Programs/ClermonTyping/bin/clermont.py", "-x", "{}/{}.xml".format(outputdir,nom)],stdout = filout)
