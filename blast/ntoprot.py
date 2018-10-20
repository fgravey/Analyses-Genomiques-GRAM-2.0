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
import sens_traduction as way
import os

def reversecomplement(seq):
    # Make reverse complement strand
    trans = str.maketrans("ATGC", "TACG")
    return seq.translate(trans)[::-1]

def strain_trad(liste,blastdir):
    print("########################################################")
    print("################# Traduction Started ##########################")
    print("########################################################")
    print("\n")
    """ DNA Traduction in protein using the Biopython library
        Inputs : list which contains all the stains you wanted working on
                 absolute path to the directory of blast
        Outputs: fasta file containing the protein sequence namde as
                 gene_nameofthestrain_protein_sequence.fasta in a dedicated
                 directory named proteins_sequences"""

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
            if fichier.find('tronquee') == -1: #not working on tronquee proteins
                ## Variables definition:
                gene = (fichier.split("/")[-1]).replace("_{}_-_nt_sequence.fasta".format(nom),"")
                prot_fasta = "{}/{}_{}_protein_sequence.fasta".\
                format(output_dir_prot,gene,nom)
                sequence_nt = []
                sequence = []

                ## Informations
                print("----> {} {}".format(nom,gene))

                ## regex definition
                regex = re.compile("^>{}_{}_".format(nom,gene))
                regex_gaps = re.compile("-")

                with open(fichier,"r") as filin:
                    for ligne in filin:
                        sequence.append(ligne[:-1])

                ### DNA Traduction regarding the way of the traduction of the gene
                ### according to the sens_traduction module

                with open(prot_fasta, "w") as filout:
                    for header in sequence:
                        if regex.search(header) and way.trad[gene] == "reverse":
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

                        elif regex.search(header) and way.trad[gene] == "strand":
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
    nom = inputdir.split("/")[-2]
    gene = query.replace("{}/".format(inputdir), "").\
    replace("_{}_protein_sequence.fasta".format(nom),"")

    subprocess.run(["blastp", "-query", "{}".format(query), "-subject", \
    "{}".format(subject), "-out", "{}{}_{}_blastp_xml.txt".format(outputdir,\
    nom, gene),'-outfmt', '5'])
    subprocess.run(["blastp", "-query", "{}".format(query), "-subject", \
    "{}".format(subject), "-out", "{}{}_{}blastp.txt".format(outputdir,nom,gene)])

def blastp_all(liste,blast_dir,database):
    #Informations
    print("########################################################")
    print("################# Blast Started ##########################")
    print("########################################################")
    print("\n")

    travail = []
    with open(liste, 'r') as filin:
        for nom in filin:
            travail.append(nom[:-1])

    for nom in travail:
        blastp_souche(nom,blast_dir,database)

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
        if subject.find('tronquee') == -1: #not working on tronquee proteins
            paths = [] ## list which will contain the path of the query and the path of the subject

            regex = re.compile("{}".format(subject.replace("{}/".format(inputdir), "").\
            replace("_{}_protein_sequence.fasta".format(nom),"")))

            ## regex is the name of the gene which will be looked for in the protein directory and
            ## to the db directory
            if regex.search(subject):
                paths.append(subject) ## adding the query absolute path paths[0]

            for query in glob.glob("{}/*.fasta".format(db)):
                if regex.search(query):
                    paths.append(query) ## adding the subject absolute directory paths[1]

        if len(paths) == 2: ## check if the path contain two paths
            blastp(paths[0], inputdir, paths[1], outputdir)

def multifasta_prot(liste,blast_dir,database):
    """ Creation of a multifasta which contains all the traducted proteins of
    the genes find in the genomes
        Inputs : list which contains all the stains you wanted working on
                 absolute path to the directory of blast
        Outpus : multifasta which contains all the traducted proteins of
        the genes find in the genomes"""

    #Informations
    print("########################################################")
    print("################# multifasta Creation ##########################")
    print("########################################################")
    print("\n")

    ##Vairables definitions
    multifasta_prot_dir = "{}/multifasta_prot".format(blast_dir)

    ### Directory Creation
    subprocess.run(["mkdir", multifasta_prot_dir])

    ## listing all the strains of the project
    travail = []
    with open(liste, 'r') as filin:
        for nom in filin:
            travail.append(nom[:-1])


    #### Adding all the genes find by blastn into a genes_list
    genes_list = []
    for nom in travail :
        inputdir = "{}/{}/proteins_sequences".format(blast_dir,nom)
        for subject in glob.glob("{}/*.fasta".format(inputdir)):
            genes = (subject.split("/")[-1]).replace("_{}_protein_sequence.fasta".format(nom),"")
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
                with open(query, 'r') as filin:
                    for ligne in filin:
                        ref_seq.append(ligne[:-1])


        ## Informations
        print("----> {}".format(g))

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
                    if subject.find('tronquee') == -1: #not working on tronquee proteins
                        if regex_2.search(subject):
                            with open(subject,'r') as filin:
                                for ligne in filin:
                                    filout.write(ligne)

def blastp_results_strain(liste,blast_dir,threshold):
    """ Parse the .xml file in order to extract specific informations such as :
            - length of the alignment
            - number of identities
            - number of positives
            - number of gaps
            - number of all the differences

        All the informations will be summarized in a variable called resultats and
        separated by ; in order to be read in a .csv file

        Inputs : blast_dir = general directory which will contain all the subdirectories
                 name = name of the strain
                 threshold =

        Outputs : resultats = list which contains all the extracted informations
        separated by a ; in order to create a .csv file
        """

    ## listing all the strains of the project
    travail = []
    with open(liste, 'r') as filin:
        for nom in filin:
            travail.append(nom[:-1])

    # Creation of two containers
    resultats = ["souche;proteine;longeur de la cible;longeur de la proteine;pourcentage de coverage;\
    pourcentage d'identite;remarques;nombre de substitutions en AA;substitutions AA\n"]
    genes_list = []
    #For each strain in the list, parsing the blastp results
    for nom in travail :
        #Variables definition:
        input_blastp_dir = "{}{}/{}_blast_protein/".format(blast_dir,nom,nom)

        #### Adding all the genes find by blastp into a genes_list
        for subject in glob.glob("{}/*_blastp_xml.txt".format(input_blastp_dir)):
            genes = (subject.split("/")[-1]).replace("{}_".format(nom),"")
            genes = genes.replace("_blastp_xml.txt","")

            if genes in genes_list:
                continue
            else :
                genes_list.append(genes)

        for gene in genes_list :
            if os.path.isfile("{}/{}_{}_blastp_xml.txt".format(input_blastp_dir,\
            nom,gene)) == True:
                #Reading the xml format of the blast results file
                with open("{}{}_{}_blastp_xml.txt".format(input_blastp_dir,\
                nom,gene)) as result_handle:
                    blast_records = NCBIXML.parse(result_handle)
                    E_VALUE_THRESH = 0.0000001
                    res = []
                    for blast_record in blast_records:
                        if blast_record.alignments:
                            for alignment in blast_record.alignments:
                                for hsp in alignment.hsps:
                                    #Variables definition
                                    souche = nom
                                    contig = blast_record.query
                                    g = (alignment.title.split(" ")[1].split("-")[0])\
                                    .split("_")[0]
                                    #query information
                                    query_start = hsp.query_start
                                    query_end = hsp.query_end

                                    #Subject information
                                    subject_longueur = int(alignment.length)
                                    sbjct_start = hsp.sbjct_start
                                    sbjct_end = hsp.sbjct_end

                                    #alignment information
                                    gaps = hsp.gaps
                                    proba = hsp.expect
                                    id = hsp.identities
                                    positive = hsp.positives
                                    alignement_longueur = len(str(hsp.query))

                                    #calculs
                                    perc_ident = int(id) / float(subject_longueur) * 100
                                    coverage = ((int(alignement_longueur) - int(gaps)) \
                                    / float(subject_longueur))
                                    perc_coverage = (((int(alignement_longueur) - int(gaps))\
                                    / float(subject_longueur)) * 100)

                                    #containers
                                    remarque = [] ## list which will contain all the 3' and/or 5' deletions
                                    substitutions = [] ## list which will contain all the nt substitutions

                                    ##keeping only sequence for which coverage is > to a threshold
                                    #by default threshold = 80
                                    if perc_coverage >= int(threshold):

                                        ###looking for N-terminal or C-terminal deletions
                                        if sbjct_start != 1 and sbjct_end != 1:
                                            remarque.append("tronquee en N-terminal")

                                        if sbjct_start != subject_longueur and sbjct_end != subject_longueur:
                                            remarque.append("tronquee en C-terminal")
                                        if not remarque :
                                            remarque.append("-")

                                        ###looking for AA substitutions
                                        for i in range(0,len(hsp.sbjct)):
                                            if hsp.sbjct[i] != hsp.query[i]:
                                                substitutions.append("{} remplace {} en position {}".format(hsp.query[i],hsp.sbjct[i], i))

                                        nb_substitutions = len(substitutions)
                                        if not substitutions:
                                            substitutions.append('-')

                                        resultats.append("{};{};{};{};{:.2f};{:.2f};{};{};{}\n"\
                                        .format(nom,gene,subject_longueur, \
                                        alignement_longueur, perc_coverage, perc_ident,\
                                        "_".join(remarque),nb_substitutions,",".join(substitutions)))

    return resultats

def blastp_results_gene(liste,blast_dir,threshold):
    """ Parse the .xml file in order to extract specific informations such as :
            - length of the alignment
            - number of identities
            - number of positives
            - number of gaps
            - number of all the differences

        All the informations will be summarized in a variable called resultats and
        separated by ; in order to be read in a .csv file

        Inputs : blast_dir = general directory which will contain all the subdirectories
                 name = name of the strain
                 threshold =

        Outputs : resultats = list which contains all the extracted informations
        separated by a ; in order to create a .csv file
        """

    ## listing all the strains of the project
    travail = []
    with open(liste, 'r') as filin:
        for nom in filin:
            travail.append(nom[:-1])

    # Creation of two containers
    genes_list = []
    dico = {}
    #For each strain in the list, parsing the blastp results
    for nom in travail :
        #Variables definition:
        input_blastp_dir = "{}{}/{}_blast_protein/".format(blast_dir,nom,nom)

        #### Adding all the genes find by blastp into a genes_list
        for subject in glob.glob("{}/*_blastp_xml.txt".format(input_blastp_dir)):
            genes = (subject.split("/")[-1]).replace("{}_".format(nom),"")
            genes = genes.replace("_blastp_xml.txt","")

            if genes in genes_list:
                continue
            else :
                genes_list.append(genes)

        for gene in genes_list :
            dico[gene] = ["souche;proteine;longeur de la cible;longeur de la proteine;pourcentage de coverage; pourcentage d'identite;remarques;nombre de substitutions en AA;substitutions AA\n"]

    #For each strain in the list, parsing the blastp results
    for nom in travail :
        #Variables definition:
        input_blastp_dir = "{}{}/{}_blast_protein/".format(blast_dir,nom,nom)

        for gene in genes_list :
            if os.path.isfile("{}/{}_{}_blastp_xml.txt".format(input_blastp_dir,\
            nom,gene)) == True:
                #Reading the xml format of the blast results file
                with open("{}/{}_{}_blastp_xml.txt".format(input_blastp_dir,\
                nom,gene)) as result_handle:
                    blast_records = NCBIXML.parse(result_handle)
                    E_VALUE_THRESH = 0.0000001
                    res = []
                    for blast_record in blast_records:
                        if blast_record.alignments:
                            for alignment in blast_record.alignments:
                                for hsp in alignment.hsps:
                                    #query information
                                    query_start = hsp.query_start
                                    query_end = hsp.query_end

                                    #Subject information
                                    subject_longueur = int(alignment.length)
                                    sbjct_start = hsp.sbjct_start
                                    sbjct_end = hsp.sbjct_end

                                    #alignment information
                                    gaps = hsp.gaps
                                    proba = hsp.expect
                                    id = hsp.identities
                                    positive = hsp.positives
                                    alignement_longueur = len(str(hsp.query))

                                    #calculs
                                    perc_ident = int(id) / float(subject_longueur) * 100
                                    coverage = ((int(alignement_longueur) - int(gaps)) \
                                    / float(subject_longueur))
                                    perc_coverage = (((int(alignement_longueur) - int(gaps))\
                                    / float(subject_longueur)) * 100)

                                    #containers
                                    remarque = [] ## list which will contain all the 3' and/or 5' deletions
                                    substitutions = [] ## list which will contain all the nt substitutions

                                    ##keeping only sequence for which coverage is > to a threshold
                                    #by default threshold = 80
                                    if perc_coverage >= int(threshold):
                                            ###looking for N-terminal or C-terminal deletions
                                        if sbjct_start != 1 and sbjct_end != 1:
                                            remarque.append("tronquee en N-terminal")

                                        if sbjct_start != subject_longueur and sbjct_end != subject_longueur:
                                            remarque.append("tronquee en C-terminal")
                                        if not remarque :
                                            remarque.append("-")

                                            ###looking for AA substitutions
                                        for i in range(0,len(hsp.sbjct)):
                                            if hsp.sbjct[i] != hsp.query[i]:
                                                substitutions.append("{} remplace {} en position {}".format(hsp.query[i],hsp.sbjct[i], i))

                                        nb_substitutions = len(substitutions)
                                        if not substitutions:
                                            substitutions.append('-')

                                        res.append("{};{};{};{};{:.2f};{:.2f};{};{};{}\n"\
                                        .format(nom,gene,subject_longueur, \
                                        alignement_longueur, perc_coverage, perc_ident,\
                                        "_".join(remarque),nb_substitutions,",".join(substitutions)))

                                        ajout = dico[gene] + res
                                        dico[gene] = ajout


    return dico

if __name__ == "__main__":

    # PARSE COMMAND LINE OPTIONS
    ##########################################################################
    parser = ArgumentParser()
    parser.add_argument("-l", "--list", dest="list", \
    help="list which contains all the name of the strains", default='')
    parser.add_argument("-b", "--blast_dir", dest="blast_directory",\
    help="Path to the directory which contains all blastn results", default='')
    parser.add_argument("-db", "--database", dest="database",\
    help="Indicate directory which contains all the fasta files as subject", default='')
    parser.add_argument("-filename", "--filename", dest="filename",\
    help="summary output file name", default='summary_blastp')
    parser.add_argument("-t", "--threshold", dest="threshold",\
    help="minimum percentage of coverage int", default='80')
    args = parser.parse_args()

    # Variables difinition
    liste = args.list
    blast_dir = args.blast_directory
    nom_fichier = args.filename
    database = args.database
    threshold = args.threshold



    # launching functions
    #strain_trad(liste,blast_dir)
    #multifasta_prot(liste,blast_dir,database)
    #blastp_all(liste,blast_dir, database)

    with open("{}{}_per_strain.csv".format(blast_dir,nom_fichier), "w") as filout:
        for ligne in blastp_results_strain(liste,blast_dir,threshold):
            filout.write(str(ligne))
        filout.write("\n")


    with open("{}{}_per_gene.csv".format(blast_dir,nom_fichier), "w") as filout:
        for clef in blastp_results_gene(liste,blast_dir,threshold).keys():
            filout.write(clef+ '\n')
            for ligne in blastp_results_gene(liste,blast_dir,threshold)[clef]:
                if len(blastp_results_gene(liste,blast_dir,threshold)[clef]) == 1:
                    filout.write("ATTENTION Proteine tronquee +++++\n")
                else:
                    filout.write(str(ligne))
            filout.write("\n")
