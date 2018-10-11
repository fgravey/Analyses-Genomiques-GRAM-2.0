#!/usr/bin/env python3
## Septembre 2018
### Gravey François

### module Loading
import subprocess
import re
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from argparse import ArgumentParser
import glob
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def reversecomplement(seq):
    # Make reverse complement strand
    trans = str.maketrans("ATGC", "TACG")
    return seq.translate(trans)[::-1]

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

def strain_trad(liste,blastdir):

    travail = []
    with open(liste, 'r') as filin:
        for nom in filin:
            travail.append(nom[:-1])

    print("Voici la liste des fichiers sur lesquels vous allez travailler", travail)
    print("Le nombre de fichier est de : {}".format(len(travail)))
    print('\n')

    for nom in travail:
        ## Variables definitions
        inputdir = "{}{}/{}_genes_sequences_fasta".format(blast_dir,nom, nom)
        output_dir_prot = "{}{}/proteins_sequences".format(blast_dir,nom)

        ## directory Creation
        subprocess.run(["mkdir", output_dir_prot])

        for fichier in glob.glob("{}/*.fasta".format(inputdir)):
            ## Variables definition:
            gene = (fichier.split("/")[-1]).split("_")[0]
            prot_fasta = "{}/{}_{}_protein_sequence.fasta".\
            format(output_dir_prot,gene,nom)
            sequence_nt = []
            sequence = []

            ## Informations
            print("########################################################")
            print("################# Working on {}-{} ##########################".format(nom,gene))
            print("########################################################")
            print("\n")

            ### Dictionary defintion
            trad = {"acrA" : "reverse", "acrB" : "reverse", "acrR" : "strand", \
            "ampC" : "reverse", "ampC2" : "strand", "ampD" : "strain",\
            "ampE" : "strand", "ampG" : "reverse", "ampH" : "reverse",\
            "ampR" : "strand", "dacA" : "strand", "dacB" : "strand",\
            "dacC" : "strand", "dacD" : "reverse", "fosA2" : "strand", \
            "lysR" : "reverse", "nagZ" : "strand", "ompC=ompK36" : "reverse",\
            "ompK35=ompF" : "reverse", "ompR" : "strand"}

            ## Regex definition
            regex_gaps = re.compile("-")

            with open(fichier,"r") as filin:
                for ligne in filin:
                    sequence.append(ligne[:-1])

            with open(prot_fasta, "w") as filout:
                for header in sequence:
                    if trad[gene] == "reverse":
                        sequence_nt = sequence[1:]
                        sequence_nt = "".join(sequence_nt)
                        if "-" in sequence_nt:
                            filout.write(">{}_{}_protein_sequence_ATTENTION_GAPS\n".format(gene,nom))
                            sequence_nt = sequence_nt.replace("-", "")
                        else:
                            filout.write(">{}_{}_protein_sequence\n".format(gene,nom))

                        coding_dna = Seq(reversecomplement(sequence_nt), generic_dna)
                        protein = coding_dna.translate(table=11, to_stop=True)
                        protein = str(protein)
                        for aa in range(0,len(protein),80):
                            filout.write(protein[aa:aa+80] + '\n')

                    elif trad[gene] == "strand":
                        sequence_nt = sequence[1:]
                        sequence_nt = "".join(sequence_nt)

                        if "-" in sequence_nt:
                            filout.write(">{}_{}_protein_sequence_ATTENTION_GAPS\n".format(gene,nom))
                            sequence_nt = sequence_nt.replace("-", "")
                        else:
                            filout.write(">{}_{}_protein_sequence\n".format(gene,nom))

                        coding_dna = Seq(sequence_nt, generic_dna)
                        protein = coding_dna.translate(table=11, to_stop=True)
                        protein = str(protein)
                        for aa in range(0,len(protein),80):
                            filout.write(protein[aa:aa+80] + '\n')

def blastp(query, inputdir, subject, outputdir):
    """Make a directory specific to each strain which contains two blast files
    first one is a .txt file and the second one is an .xml file.
    The xml file will be parse by bipython
    Inputs : nom = name of the strain
             outputdir = directory which will contains all the subdirectories

    Outputs : a specific output by strain which contains the two fasta files
    """
    #Variable definition
    nom = query.replace("{}/".format(inputdir), "").split("_")[1]
    gene = query.replace("{}/".format(inputdir), "").split("_")[0]

    subprocess.run(["blastp", "-query", "{}".format(query), "-subject", \
    "{}".format(subject), "-out", "{}/{}_{}_blastp_xml.txt".format(outputdir,\
    nom, gene),'-outfmt', '5'])
    subprocess.run(["blastp", "-query", "{}".format(query), "-subject", \
    "{}".format(subject), "-out", "{}/{}_{}blastp.txt".format(outputdir,nom,gene)])

def blastp_souche(nom, inputdir, db):
    """ function which blastp all the proteins created against the reference proteins
    using the blastp function
        Input : nom = name of the strains
                inputdir = path to the general blast project
                db = path to the directory which contains all the reference proteins

        Outputs : two blastp results file one in xml format and the second one
        in .txt format. All the files are located in the nom_blast_protein"""

    ## Variables definition and directory creation
    outputdir = "{}/{}/{}_blast_protein".format(inputdir,nom,nom)
    subprocess.run(["mkdir", outputdir])
    inputdir = "{}/{}/proteins_sequences".format(inputdir,nom)

    ## Looking for all the fasta files present in the protein directory
    for subject in glob.glob("{}/*.fasta".format(inputdir)):
        paths = [] ## list which will contain the path of the query and the path of the subject

        regex = re.compile("{}_".format(subject.replace("{}/".format(inputdir), "").split("_")[0]))
        ## regex is the name of the gene which will be looked for in the protein directory and
        ## to the db directory
        if regex.search(subject):
            print(subject)
            paths.append(subject) ## adding the query absolute path paths[0]

        for query in glob.glob("{}/*.fasta".format(db)):
            if regex.search(query):
                print(query)
                paths.append(query) ## adding the subject absolute directory paths[1]

        print(paths)
        print(len(paths))

        if len(paths) == 2: ## check if the path contain two paths
            blastp(paths[0], inputdir, paths[1], outputdir)

if __name__ == "__main__":

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
    help="Indicate directory which contains all the fasta files as subject", default='')
    parser.add_argument("-filename", "--filename", dest="filename",\
    help="summary output file name", default='summary_blast')
    parser.add_argument("-t", "--threshold", dest="threshold",\
    help="minimum percentage of coverage int", default='80')
    args = parser.parse_args()

    # Variables difinition
    liste = args.list
    blast_dir = args.blast_directory
    outputdir = args.out_path
    fasta_extension = args.extension
    nom_fichier = args.filename
    database = args.database
    threshold = args.threshold

    strain_trad(liste,blast_dir)
    #blastp_souche("FAY1", blast_dir, database)
    #liste,blastdir,database

    ##Vairables definitions
    multifasta_prot_dir = "{}/multifasta_prot".format(blast_dir)

    ### Directory Creation
    subprocess.run(["mkdir", multifasta_prot_dir])

    ## listing all the strains of the project
    travail = []
    with open(liste, 'r') as filin:
        for nom in filin:
            travail.append(nom[:-1])

    print("Voici la liste des fichiers sur lesquels vous allez travailler", travail)
    print("Le nombre de fichier est de : {}".format(len(travail)))
    print('\n')

    #### Adding all the genes find by blastn into a genes_list
    genes_list = []
    for nom in travail :
        inputdir = "{}/{}/proteins_sequences".format(blast_dir,nom)
        for subject in glob.glob("{}/*.fasta".format(inputdir)):
            genes = (subject.split("/")[-1]).split("_")[0]
            if genes in genes_list:
                continue
            else :
                genes_list.append(genes)

    ### parsing the files in order to extract the proteins sequences
    for g in genes_list:
        ### Regex definition:
        regex = re.compile("{}_".format(g))


        ## Variable definition
        ref_seq = [] # list which will contain the header and the sequence of the reference protein
        for query in glob.glob("{}/*.fasta".format(database)):
            if regex.search(query):
                print(query)
                with open(query, 'r') as filin:
                    for ligne in filin:
                        ref_seq.append(ligne[:-1])

        ## Writing the reference sequence into the multifasta file
        with open("{}/{}_multifasta_protein.fasta".format(multifasta_prot_dir,g),'w')\
        as filout:
            filout.write(ref_seq[0]+ "_reference_sequence" + '\n')
            ref_seq = "".join(ref_seq[1:])
            for i in range(0,len(ref_seq),80):
                filout.write(ref_seq[i:i+80]+'\n')

        with open("{}/{}_multifasta_protein.fasta".format(multifasta_prot_dir,g),'a')\
        as filout:
            for nom in travail:
                ## Regex definition
                regex_2 = re.compile("{}_{}".format(g,nom))
                ## Variable definition
                inputdir = "{}/{}/proteins_sequences".format(blast_dir,nom)
                ## For each strain looking for the protein sequence
                for subject in glob.glob("{}/*.fasta".format(inputdir)):
                    if regex_2.search(subject):
                        print(subject)
                        with open(subject,'r') as filin:
                            for ligne in filin:
                                filout.write(ligne)




    # for query in glob.glob("{}/*.fasta".format(database)):
        # gene = (query.split("/")[-1]).split("_")[0]
