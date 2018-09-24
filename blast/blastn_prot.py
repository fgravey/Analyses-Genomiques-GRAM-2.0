#!/usr/bin/env python3
## Septembre 2018
### Gravey François

### module Loading
import subprocess
import re
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
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

def gene_to_protein(nom,inputdir):
    """Parse the _blast_xml result file and extract all the nucleotidique
    sequences find. Then the nucleotidique sequence is translated using the
    traduction function and write in a specific file
        Input : nom : name of the strain
                inputdir : general directory of the blast project

        Outputs : fasta files containing the proteic sequence of all the genes
        find in the query genome"""

    #Variables definition
    inputdir = "{}/{}".format(inputdir, nom)
    outputdir = "{}/protein".format(inputdir)
    subprocess.run(["mkdir", outputdir])

    #Parsing the _blast_xml.txt file
    with open("{}/{}_blast_xml.txt".format(inputdir,nom), "r") as result_handle:
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
                        a_traduire = sequence.replace("-","")
                        protein = "".join(traduction(a_traduire.lower())).upper()
                        with open("{}/{}_{}_{}_protein.fasta".format(outputdir,\
                        gene,nom,contig), "w") as filout:
                            filout.write(">{}-{}-{}-protein\n".format(souche, \
                            contig, gene))
                            for i in range(0,len(protein),80):
                                filout.write(protein[i:i+80]+'\n')

def blastp(query, inputdir, subject, outputdir):
    """Make a directory specific to each strain which contains two blast files
    first one is a .txt file and the second one is an .xml file.
    The xml file will be parse by bipython
    Inputs : nom = name of the strain
             outputdir = directory which will contains all the subdirectories

    Outputs : a specific output by strain which contains the two fasta files
    """
    nom = query.replace("{}/".format(inputdir), "").split("_")[0]
    contig = query.replace("{}/".format(inputdir), "").split("_")[5]

    subprocess.run(["blastp", "-query", "{}".format(query), "-subject", \
    "{}".format(subject), "-out", "{}/{}_contig_{}_blastp_xml.txt".format(outputdir,\
    nom, contig),'-outfmt', '5'])
    subprocess.run(["blastp", "-query", "{}".format(query), "-subject", \
    "{}".format(subject), "-out", "{}/{}_contig_{}blastp.txt".format(outputdir,nom,contig)])


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
    inputdir = "{}/{}/protein".format(inputdir,nom)

    ## Looking for all the fasta files present in the protein directory
    for fichier in glob.glob("{}/*.fasta".format(inputdir)):
        paths = [] ## list which will contain the path of the query and the path of the subject

        regex = re.compile("{}_".format(fichier.replace("{}/".format(inputdir), "").split("_")[0]))
        ## regex is the name of the gene which will be looked for in the protein directory and
        ## to the db directory

        paths.append(fichier) ## adding the query absolute path paths[0]

        for fichier in glob.glob("{}/*.fasta".format(db)):
            if regex.search(fichier):
                paths.append(fichier) ## adding the subject absolute directory paths[1]

        if len(paths) == 2: ## check if the path contain two paths
            blastp(paths[0], inputdir, paths[1], outputdir)

#inputdir = "/Users/Francois/Desktop/essai_blast"
#nom = "FAY1"
#db = "/Users/Francois/blast_data_base/ecloacae/proteines_cluster_6"
#blastp_souche(nom, inputdir, db)

for fichier in glob.glob("{}/*_xml.txt".format("/Users/Francois/Desktop/essai_blast/FAY1/FAY1_blast_protein")):
    with open(fichier, "r") as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        E_VALUE_THRESH = 0.04
        for blast_record in blast_records:
            if blast_record.alignments:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        ## Variables definition
                        #souche = nom
                        contig = blast_record.query
                        gene = alignment.title.split(" ")[1].split("-")[0]
                        diff = int(alignment.length)-(int(hsp.identities)+int(hsp.gaps))
                        longueur = alignment.length
                        proba = hsp.expect
                        id = hsp.identities
                        positive = hsp.positives
                        gaps = hsp.gaps
                        sequence = hsp.query
                        match = hsp.match
                        subject = hsp.sbjct
                        print(contig)
                        print(gene)
                        #print(diff)
                        #print(sequence)
                        #print(subject)
                        print(match)
                        #for i in range(0,len(subject)):
                            #if sequence[i] != subject[i]:
                                #print("{} remplace {} en position {}".format(sequence[i], subject[i], i))
            else:
                contig = blast_record.query
                print("{} no match !!".format(contig))





#gene_to_protein('FAY1','/Users/Francois/Desktop/essai_blast')
