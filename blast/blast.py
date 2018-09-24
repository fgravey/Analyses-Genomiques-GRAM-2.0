#!/usr/bin/env python3
## Aout 2018
### Gravey François

### module Loading
import subprocess
import re
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from argparse import ArgumentParser

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
    outputdir = "{}{}".format(outputdir,nom)
    subprocess.run(["mkdir", outputdir])
    #launch blast between all the files in the list and the database
    fasta = "{}{}{}".format(fasta_dir,nom,fasta_extension)
    subprocess.run(["blastn", "-query", "{}".format(fasta), "-db", \
    "{}".format(database), "-out", "{}{}_blast_xml.txt".format(outputdir,nom),\
    '-outfmt', '5'])
    subprocess.run(["blastn", "-query", "{}".format(fasta), "-db", \
    "{}".format(database), "-out", "{}{}_blast.txt".format(outputdir,nom)])


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
        separated by a ;
        """
    ## Variables definition
    inputdir = "{}{}".format(outputdir,nom)
    multifasta_nt_sequence_dir = "{}/{}".format(outputdir,nom)
    with open("{}/{}_blast_xml.txt".format(outputdir,nom)) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        E_VALUE_THRESH = 0.0000001
        compteur = 0 #count the number of genes found in the strain's genome
        resultats = []
        for blast_record in blast_records:
            #print(dir(blast_record))
            if blast_record.alignments:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        ## Variables definition
                        ## General information
                        souche = nom
                        contig = blast_record.query
                        gene = alignment.title.split(" ")[1].split("-")[0]
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
                        compteur = compteur + 1
                        perc_ident = int(id) / float(subject_longueur) * 100
                        coverage = ((int(alignement_longueur) - int(gaps))
                                    / float(subject_longueur))
                        perc_coverage = (((int(alignement_longueur) - int(gaps))
                                          / float(subject_longueur)) * 100)

                        ### keeping only sequence for which coverage is > to a threshold
                        ## by default threshold = 80
                        if perc_coverage >= "{}".format(threshold):
                            print(gene)
                            print(subject_longueur)
                            print((int(query_end)+1)-int(query_start))
                            print("coverage")
                            print(perc_coverage)
                            print("identity")
                            print(perc_ident)
                            print(sbjct_start, sbjct_end)

                            ###looking for 3' or 5' deletions
                            if sbjct_start != 1 and sbjct_end != 1:
                                print("tronquee en 5'")
                            if sbjct_start != subject_longueur and sbjct_end != subject_longueur:
                                print("tronquee en 3'")

                            ###looking for which DNA strand the gene is located
                            if sbjct_start > sbjct_end:
                                print(">{}_sequence".format(nom))
                                for i in range(0,len(hsp.query),80):
                                    print(reversecomplement(hsp.query)[i:i+80])

                                print(">Ref_qeuence")
                                for i in range(0,len(hsp.sbjct),80):
                                    print(reversecomplement(hsp.sbjct)[i:i+80])
                                for i in range(0,len(reversecomplement(hsp.sbjct))):
                                    if reversecomplement(hsp.sbjct)[i] != reversecomplement(hsp.query)[i]:
                                        print("{} remplace {} en position {}".format(reversecomplement(hsp.query)[i],reversecomplement(hsp.sbjct)[i], i))
                            else:
                                print(">{}_sequence".format(nom))
                                print(hsp.query)
                                print("######################################")
                                print(">Ref_qeuence")
                                print(hsp.sbjct)
                                for i in range(0,len(hsp.sbjct)):
                                    if hsp.sbjct[i] != hsp.query[i]:
                                        print("{} remplace {} en position {}".format(reversecomplement(hsp.query)[i],reversecomplement(hsp.sbjct)[i], i))
                            print("\n")

                        resultats.append("{};{};{};{};{};{};{}\n".format(souche,\
                        gene, contig,subject_longueur, id, positive, gaps))
    resultats.append("le nombre de genes trouves est de {}\n".format(compteur))
    resultats = "".join(resultats)
    return resultats

def blast_nt_result_filout(liste, fasta_dir, outputdir, fasta_extension, database):
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
    sortie = ["souche;gene;contig;longeur;nb_identites;nb_positives;nb_gaps;nb_differences\n"]
    for nom in travail:
        blastn(nom,outputdir,fasta_dir)
        sortie.append(blast_nt_result(nom,outputdir))

    return sortie


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
# Variables difinition
liste = args.list
fasta_dir = args.fasta_path_dir
outputdir = args.out_path
fasta_extension = args.extension
nom_fichier = args.filename
database = args.database
threshold = args.threshold

#with open("{}/{}.csv".format(outputdir,nom_fichier), 'w') as filout:
    #for results in blast_nt_result_filout(liste, fasta_dir, outputdir, fasta_extension, database):
        #filout.write(results)


fasta_extension = ".agp.fasta"
database = "/Users/Francois/blast_data_base/ecloacae/interet/genes_ecc_fg.fasta"
blast_nt_result("FAY1", "/Users/Francois/Desktop/essai_blast")
#blastn("FAY1", "/Users/Francois/Desktop/essai_blast", "/Users/Francois/Desktop")
