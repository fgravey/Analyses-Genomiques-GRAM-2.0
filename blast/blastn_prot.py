#!/usr/bin/env python3
## Aout 2018
### Gravey François

### module Loading
import subprocess
import re
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML

##Function defintion
def traduction(liste):
	"""Fonction qui traduit une séquence d'ADN issue d'un gène Fasta en une
	séquence protéique"""

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

def gene_to_protein(nom,outputdir):
    outputdir = "{}/{}".format(outputdir,nom)
    with open("{}/{}_blast_xml.txt".format(outputdir,nom)) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        E_VALUE_THRESH = 0.04
        for blast_record in blast_records:
            if blast_record.alignments:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        ## Variables definition
                        souche = nom
                        contig = blast_record.query
                        gene = alignment.title.split(" ")[1].split("-")[0]
                        diff = int(alignment.length)-(int(hsp.identities)+int(hsp.gaps))
                        longueur = alignment.length
                        proba = hsp.expect
                        id = hsp.identities
                        positive = hsp.positives
                        gaps = hsp.gaps
                        sequence = hsp.query
                        protein = "".join(traduction(sequence.lower())).upper()
                        with open("{}/{}_{}_protein.fasta".format(outputdir,\
                        gene,nom), "w") as filout:
                            filout.write(">{}-{}-{}-protein\n".format(souche, \
                            contig, gene))
                            for i in range(0,len(protein),80):
                                filout.write(protein[i:i+80]+'\n')

def blastp(nom,gene,outputdir,database):
    """Make a directory specific to each strain which contains two blast files
    first one is a .txt file and the second one is an .xml file.
    The xml file will be parse by bipython
    Inputs : nom = name of the strain
             outputdir = directory which will contains all the subdirectories

    Outputs : a specific output by strain which contains the two fasta files
    """
    outputdir = "{}/{}".format(outputdir,nom)
    fasta = "{}/{}_{}_protein.fasta".format(outputdir,gene,nom)
    database = "{}/{}.fasta"

    #launch blast between all the files in the list and the database
    subprocess.run(["blastp", "-query", "{}".format(fasta), "-db", \
    "{}".format(database), "-out", "{}/{}_{}_blastp_xml.txt".format(outputdir,\
    gene,nom),'-outfmt', '5'])
    subprocess.run(["blastp", "-query", "{}".format(fasta), "-db", \
    "{}".format(database), "-out", "{}/{}_{}_blastp.txt".format(outputdir,gene,nom)])

#gene_to_protein('201601729','/Users/Francois/Desktop/essai_blast')
