#!/usr/bin/env python3
## Aout 2018
### Gravey François

### module Loading
import subprocess
import re
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from argparse import ArgumentParser
import glob

##Function defintion
def reversecomplement(seq):
    # Make reverse complement strand
    trans = str.maketrans("ATGC", "TACG")
    return seq.translate(trans)[::-1]

def blastn(nom,outputdir,fasta_dir):
    """Make a directory specific to each strain wich contains two blast files
    first one is a .txt file and the second one is an .xml file.
    The xml file will be parse by bipython
    Inputs : nom = name of the strain
             outputdir = directory which will contains all the subdirectories
             fasta_dir = directory which contains the fasta files

    Outputs : a specific output by strain which contains the two fasta files
    """
    outputdir_blast = "{}/{}".format(outputdir,nom)
    subprocess.run(["mkdir", outputdir_blast])
    #launch blast between all the files in the list and the database
    fasta = "{}{}{}".format(fasta_dir,nom,fasta_extension)
    subprocess.run(["blastn", "-query", "{}".format(fasta), "-db", \
    "{}".format(database), "-out", "{}/{}_blast_xml.txt".format(outputdir_blast,nom),\
    '-outfmt', '5'])
    subprocess.run(["blastn", "-query", "{}".format(fasta), "-db", \
    "{}".format(database), "-out", "{}/{}_blast.txt".format(outputdir_blast,nom)])


def blast_nt_result(nom,inputdir, threshold):
    """ Parse the .xml file in order to extract specific informations such as :
            - contig in which the gene of interest have been recovered
            - name of the gene of interest
            - length of the alignment
            - number of identities
            - number of positives
            - number of gaps
            - number of all the differences

        All the informations will be summarized in a variable called resultats and
        separated by ; in order to be read in a .csv file

        Inputs : outputdir = general directory which will contain all the subdirectories
                 name = name of the strain

        Outputs : resultats = list which contains all the extracted informations
        separated by a ; in order to create a .csv file
        """
    ## Variables definition
    inputdir = "{}{}".format(inputdir,nom)
    multifasta_nt_sequence_dir = "{}/{}_genes_sequences_fasta".format(inputdir,nom)
    subprocess.run(["mkdir", multifasta_nt_sequence_dir])

    ## Reading the xml format of the blast results file
    with open("{}/{}_blast_xml.txt".format(inputdir,nom)) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        E_VALUE_THRESH = 0.0000001
        resultats = []
        compteur = 0 #count the number of genes found in the strain's genome
        for blast_record in blast_records:
            #print(dir(blast_record))
            if blast_record.alignments:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        ## Variables definition
                        ## General information
                        souche = nom
                        contig = blast_record.query
                        gene = (alignment.title.split(" ")[1].split("-")[0])

                        ## query information
                        query_start = hsp.query_start
                        query_end = hsp.query_end

                        ## Subject information
                        subject_longueur = int(alignment.length)
                        sbjct_start = hsp.sbjct_start
                        sbjct_end = hsp.sbjct_end

                        ## alignment information
                        gaps = hsp.gaps
                        proba = hsp.expect
                        id = hsp.identities
                        positive = hsp.positives
                        alignement_longueur = len(str(hsp.query))

                        ## calculs
                        perc_ident = int(id) / float(subject_longueur) * 100
                        coverage = ((int(alignement_longueur) - int(gaps))
                                    / float(subject_longueur))
                        perc_coverage = (((int(alignement_longueur) - int(gaps))
                                          / float(subject_longueur)) * 100)

                        ## containers
                        remarque = [] ## list which will contain all the 3' and/or 5' deletions
                        substitutions = [] ## list which will contain all the nt substitutions

                        ### keeping only sequence for which coverage is > to a threshold
                        ## by default threshold = 80
                        if perc_coverage >= int(threshold):
                            compteur = compteur + 1

                            ###looking for 3' or 5' deletions
                            if sbjct_start != 1 and sbjct_end != 1:
                                remarque.append("tronquee en 5'")
                            if sbjct_start != subject_longueur and sbjct_end != subject_longueur:
                                remarque.append("tronquee en 3'")
                            if not remarque :
                                remarque.append("-")

                            with open("{}/{}_{}_{}_nt_sequence.fasta".format(multifasta_nt_sequence_dir,gene,nom,remarque[0].split(" ")[0]), "w") as filout:
                                ###looking for which DNA strand the gene is located
                                if sbjct_start > sbjct_end:
                                    dna_strand = "-1"
                                    filout.write(">{}_{}_sequence_reverse".format(nom,gene)+'\n')
                                    for i in range(0,len(hsp.query),80):
                                        filout.write(reversecomplement(hsp.query)[i:i+80]+'\n')

                                    #filout.write(">{}_Ref_qeuence".format(gene)+'\n')
                                    #for i in range(0,len(hsp.sbjct),80):
                                        #filout.write(reversecomplement(hsp.sbjct)[i:i+80]+'\n')

                                    for i in range(0,len(reversecomplement(hsp.sbjct))):
                                        if reversecomplement(hsp.sbjct)[i] != reversecomplement(hsp.query)[i]:
                                            substitutions.append("{} remplace {} en position {}".format(reversecomplement(hsp.query)[i],reversecomplement(hsp.sbjct)[i], i))
                                else:
                                    dna_strand = "1"
                                    filout.write(">{}_{}_sequence".format(nom,gene)+'\n')
                                    for i in range(0,len(hsp.query),80):
                                        filout.write(hsp.query[i:i+80]+'\n')

                                    #filout.write(">{}_Ref_sequence".format(gene)+'\n')
                                    #for i in range(0,len(hsp.sbjct),80):
                                        #filout.write(hsp.sbjct[i:i+80]+'\n')

                                    for i in range(0,len(hsp.sbjct)):
                                        if hsp.sbjct[i] != hsp.query[i]:
                                            substitutions.append("{} remplace {} en position {}".format(hsp.query[i],hsp.sbjct[i], i))

                            nb_substitutions = len(substitutions)
                            if not substitutions:
                                substitutions.append('-')

                            resultats.append("{};{};{};{};{};{:.2f};{:.2f};{};{};{};{}\n"\
                            .format(souche,gene, contig,subject_longueur, \
                            alignement_longueur, perc_coverage, perc_ident,\
                            dna_strand,"_".join(remarque),nb_substitutions,",".join(substitutions)))

    resultats.append("{};le nombre de genes trouves est de {};{};{};{};{};{};{};{};{};{}\n".\
    format(nom, compteur, "-", "-", "-", "-", "-", "-", "-", "-", "-"))
    resultats = "".join(resultats)

    return resultats

def blast_nt_result_filout(liste, fasta_dir, outputdir, fasta_extension, database, threshold):
    """ Parse a .txt file which contains all the names of the strains and launch
    the blastn function and the blast_nt_result function

    Inputs : liste = .txt file which contains all the name of the strains
             fasta_dir = directory which contains the fasta files
             outputdir = general directory which will contain all the subdirectories
             fasta_extension = information about recursive information about the
             fasta extension such as [name].XXXXX.fasta DO NOT ADD THE FIRST POINT
             database = absolute path to the fasta file which have been the source
             of the creation of the database

    Outputs : .cvs file which contains all the informations extracted. The file
    will take place into the outputdir directory"""

    travail = []
    with open(liste, 'r') as filin:
        for nom in filin:
            travail.append(nom[:-1])

    print("Voici la liste des fichiers sur lesquels vous allez travailler", travail)
    print("Le nombre de fichier est de : {}".format(len(travail)))
    sortie = ["souche;gene;contig;longeur de la cible;longeur du gene;pourcentage\
     de coverage;pourcentage d'identite;dna_strand;remarques;nombre de substitutions nucleotidiques;substitutions nucleotidiques\n"]
    for nom in travail:
        print("############################################################################")
        print("########################### Working on {} ###############################".format(nom))
        print("############################################################################")
        print('\n')
        blastn(nom,outputdir,fasta_dir)
        sortie.append(blast_nt_result(nom,outputdir, threshold))

    return sortie

def multifasta_nt(liste, outputdir, database):
    """ Parse a .txt file which contains all the names of the strains and launch
        the blastn function and the blast_nt_result function

        Inputs : liste = .txt file which contains all the name of the strains
                 outputdir = general directory which will contain all the subdirectories

        Outputs : multi-sequence fasta file, one per gene, which contains all the sequences
    find in each strain """

    travail = [] #list which will contain all the strains name
    input = outputdir
    with open(liste, 'r') as filin:
        for nom in filin:
            travail.append(nom[:-1]) #adding to travail list all the strain names

    genes = [] #empty list which will contain all the genes name find by blast whit the first strain
    for nom in travail :
        for fichier in glob.glob("{}/{}/{}_genes_sequences_fasta/*.fasta".format(input, nom,nom)):
            if fichier.find('tronquee') == -1: #not working on tronquee proteins
                g = (fichier.split("/")[-1]).replace('_{}_-_nt_sequence.fasta'.format(nom),'')
                if g in genes:
                    continue
                else :
                    genes.append(g)

    multifasta_gene_output = "{}/multifasta_nt".format(input) ## Creation of a directory
    subprocess.run(["mkdir", multifasta_gene_output])
    print(genes)
    for g in genes:

        #Informations
        print("########################################################")
        print("################# Working on {} multifasta ##########################".format(g))
        print("########################################################")
        print("\n")

        # regex definition
        regex = re.compile("{}_".format(g)) #Name of the genes
        regex_2 = re.compile("^>") #Looking for header in fasta file
        regex_3 = re.compile("^>{}[_|-]".format(g))

        with open("{}/{}_multifasta.fasta".format(multifasta_gene_output, g), "w") as filout:
            with open("{}".format(database), "r") as filin:
                ligne = filin.readline()
                while regex_3.search(ligne) == None:
                    ligne = filin.readline()

                if regex_3.search(ligne):
                    filout.write(ligne[:-1]+'_reference_sequence\n')
                    ligne = filin.readline()

                    sequence = []
                    while regex_2.search(ligne) == None and ligne != '':
                        sequence.append(ligne[:-1])
                        ligne = filin.readline()

                    sequence = ''.join(sequence)
                    for i in range(0,len(sequence),80):
                        filout.write(sequence[i:i+80]+'\n')

        with open("{}/{}_multifasta.fasta".format(multifasta_gene_output, g), "a") as filout:
            for nom in travail:
                for fichier in glob.glob("{}/{}/{}_genes_sequences_fasta/*.fasta".format(input, nom,nom)):
                    if fichier.find('tronquee') == -1:#not working on tronquee proteins
                        if regex.search(fichier):
                            with open("{}".format(fichier), "r") as filin:
                                ligne = filin.readline()
                                if regex_2.search(ligne):
                                    filout.write(ligne)
                                    ligne = filin.readline()
                                    while regex_2.search(ligne) == None and ligne != '':
                                        filout.write(ligne)
                                        ligne = filin.readline()


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
parser.add_argument("-f", "--fasta_path", dest="fasta_path_dir",\
help="Path to the directory which contains all the fasta to analyse", default='')
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
    fasta_dir = args.fasta_path_dir
    outputdir = args.out_path
    fasta_extension = args.extension
    nom_fichier = args.filename
    database = args.database
    threshold = args.threshold

    with open("{}{}.csv".format(outputdir,nom_fichier), 'w') as filout:
        for results in blast_nt_result_filout(liste, fasta_dir, outputdir,\
         fasta_extension, database, threshold):
            filout.write(results)

    multifasta_nt(liste,outputdir,database)
