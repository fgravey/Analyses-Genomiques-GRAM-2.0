#!/usr/bin/env python3
## Septembre 2018
### Gravey François

### module Loading
import subprocess
import re
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from argparse import ArgumentParser
import glob




##Function defintion
def traduction(liste):
	"""Fonction qui traduit une séquence d'ADN issue d'un gène Fasta en une
	séquence protéique
        Input : liste qui contient la sequence nucleotidique a a_traduire
        Ouput : liste qui contient la sequence proteique traduite"""

	prot= []
	seq = liste
	dico = {}
	gencode = {
    'ata':'i', 'atc':'i', 'att':'i', 'atg':'m',
    'aca':'t', 'acc':'t', 'acg':'t', 'act':'t',
    'aac':'n', 'aat':'n', 'aaa':'k', 'aag':'k',
    'agc':'s', 'agt':'s', 'aga':'r', 'agg':'r',
    'cta':'l', 'ctc':'l', 'ctg':'l', 'ctt':'l',
    'cca':'p', 'ccc':'p', 'ccg':'p', 'cct':'p',
    'cac':'h', 'cat':'h', 'caa':'q', 'cag':'q',
    'cga':'r', 'cgc':'r', 'cgg':'r', 'cgt':'r',
    'gta':'v', 'gtc':'v', 'gtg':'v', 'gtt':'v',
    'gca':'a', 'gcc':'a', 'gcg':'a', 'gct':'a',
    'gac':'d', 'gat':'d', 'gaa':'e', 'gag':'e',
    'gga':'g', 'ggc':'g', 'ggg':'g', 'ggt':'g',
    'tca':'s', 'tcc':'s', 'tcg':'s', 'tct':'s',
    'ttc':'f', 'ttt':'f', 'tta':'l', 'ttg':'l',
    'tac':'y', 'tat':'y', 'taa':'_', 'tag':'_',
    'tgc':'c', 'tgt':'c', 'tga':'_', 'tgg':'w'}


	debut = seq.find('atg')
	for i in range(debut, len(seq)-2, 3):
		mot = seq[i:i+3]
		if mot != 'taa' and mot != 'tga' and mot != 'tag':
			prot.append(gencode[mot])
		else:
			break


	return prot


####################################################################################
################################# Main #############################################
####################################################################################
# PARSE COMMAND LINE OPTIONS
##########################################################################
parser = ArgumentParser()
parser.add_argument("-l", "--list", dest="list", \
help="list which contains all the name of the strains", default='')
parser.add_argument("-o", "--outputPath", dest="out_path",\
help="Path to blast output", default='')
parser.add_argument("-b", "--blast_dir", dest="blast_directory",\
help="Path to the directory which contains all blastn results", default='')
parser.add_argument("-e", "--extension", dest="extension",\
help="fasta extension such as .fasta .fa .awked.fasta .agp.fasta", default='')
parser.add_argument("-db", "--database", dest="database",\
help="Indicate the name of the fasta file used to create the database", default='')
parser.add_argument("-filename", "--filename", dest="filename",\
help="summary output file name", default='summary_blast')
parser.add_argument("-t", "--threshold", dest="threshold",\
help="minimum percentage of coverage int", default='80')
args = parser.parse_args()

###############################################################################

if __name__ == "__main__":
    # Variables difinition
    liste = args.list
    blast_dir = args.blast_directory
    outputdir = args.out_path
    fasta_extension = args.extension
    nom_fichier = args.filename
    database = args.database
    threshold = args.threshold

#liste,blastdir
travail = []
with open(liste, 'r') as filin:
    for nom in filin:
        travail.append(nom[:-1])

print("Voici la liste des fichiers sur lesquels vous allez travailler", travail)
print("Le nombre de fichier est de : {}".format(len(travail)))
for nom in travail:
    ## Variables definitions
    inputdir = "{}{}/{}_genes_sequences_fasta".format(blast_dir,nom, nom)
    output_dir_prot = "{}{}/proteins_sequences".format(blast_dir,nom)

    ## directory Creation
    subprocess.run(["mkdir", output_dir_prot])

    for fichier in glob.glob("{}/*.fasta".format(inputdir)):
        ## Variables definition:
        gene = (fichier.split("/")[-1]).split("_")[0]
        print(gene)
        prot_fasta = "{}/{}_{}_protein_sequence.fasta".\
        format(output_dir_prot,gene,nom)

        ## regex definition
        regex = re.compile("^>{}_{}_".format(nom,gene))

        print(prot_fasta)
        with open(fichier,"r") as filin:
            for ligne in filin:
                if regex.search(ligne):
                    header = ">{}_{}_protein_sequence".format(nom,gene)
                    print(header)
                else:
                    print(ligne)
