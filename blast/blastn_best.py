#!/usr/bin/env python3
## Novembre 2018
### Gravey François

import subprocess
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
import glob

def reversecomplement(seq):
    # Make reverse complement strand
    trans = str.maketrans("ATGC", "TACG")
    return seq.translate(trans)[::-1]

def blastn(nom,outputdir,fasta_dir,fasta_extension,database):
    """Make a directory specific to each strain wich contains two blast files
    first one is a .txt file and the second one is an .xml file.
    The xml file will be parse by bipython
    Inputs : nom = name of the strain
             outputdir = directory which will contains all the subdirectories
             fasta_dir = directory which contains the fasta files

    Outputs : a specific output by strain which contains the two fasta files
    """

    #launch blast between all the files in the list and the database
    fasta = "{}{}{}".format(fasta_dir,nom,fasta_extension)
    subprocess.run(["blastn", "-query", "{}".format(fasta), "-db", \
    "{}".format(database), "-out", "{}{}_blast_xml.txt".format(outputdir,nom),\
    '-outfmt', '5'])
    subprocess.run(["blastn", "-query", "{}".format(fasta), "-db", \
    "{}".format(database), "-out", "{}{}_blast.txt".format(outputdir,nom)])


def blast_nt_best_result(nom,inputdir, threshold):
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

    ## Reading the xml format of the blast results file
    with open("{}{}_blast_xml.txt".format(inputdir,nom)) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        E_VALUE_THRESH = 0.0000001
        best = []
        for blast_record in blast_records:
            if blast_record.alignments:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        #Variables definition
                        #General information
                        souche = nom
                        contig = blast_record.query
                        gene = (alignment.title.split(" ")[1].split("-")[0])

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
                        coverage = ((int(alignement_longueur) - int(gaps))
                                    / float(subject_longueur))
                        perc_coverage = (((int(alignement_longueur) - int(gaps))
                                          / float(subject_longueur)) * 100)

                        #containers
                        remarque = [] ## list which will contain all the 3' and/or 5' deletions
                        substitutions = [] ## list which will contain all the nt substitutions

                        ##keeping only sequence for which coverage is > to a threshold
                        #by default threshold = 80
                        if perc_coverage >= int(threshold):
                            ###looking for 3' or 5' deletions
                            if sbjct_start != 1 and sbjct_end != 1:
                                remarque.append("tronquee en 5'")
                            if sbjct_start != subject_longueur and sbjct_end != subject_longueur:
                                remarque.append("tronquee en 3'")
                            if not remarque :
                                remarque.append("-")

                            sequence = []
                            ###looking for which DNA strand the gene is located
                            if sbjct_start > sbjct_end:
                                dna_strand = "-1"
                                sequence = reversecomplement(hsp.query)
                                for i in range(0,len(reversecomplement(hsp.sbjct))):
                                    if reversecomplement(hsp.sbjct)[i] != reversecomplement(hsp.query)[i]:
                                        substitutions.append("{} remplace {} en position {}".format(reversecomplement(hsp.query)[i],reversecomplement(hsp.sbjct)[i], i))
                            else:
                                dna_strand = "1"
                                sequence = hsp.query
                                for i in range(0,len(hsp.sbjct)):
                                    if hsp.sbjct[i] != hsp.query[i]:
                                        substitutions.append("{} remplace {} en position {}".format(hsp.query[i],hsp.sbjct[i], i))

                            nb_substitutions = len(substitutions)
                            if not substitutions:
                                substitutions.append('-')

                            if not best:
                                best = "{};{};{};{};{};{:.2f};{:.2f};{:.2f};{:.2f};{};{};{};{}\n".\
                                format(nom,gene,contig,subject_longueur,\
                                alignement_longueur,perc_ident,perc_coverage,\
                                hsp.score,hsp.bits, "".join(remarque), \
                                nb_substitutions,"".join(substitutions),sequence)

                            elif (float(hsp.bits) > float(best.split(";")[8])):
                                best = "{};{};{};{};{};{:.2f};{:.2f};{:.2f};{:.2f};{};{};{};{}\n".\
                                format(nom,gene,contig,subject_longueur,\
                                alignement_longueur,perc_ident,perc_coverage,\
                                hsp.score,hsp.bits, "".join(remarque), \
                                nb_substitutions,"".join(substitutions),sequence)

    #If there is no match during the blastn process
    if not best:
        best = "{0};ND;{1};{1};{1};{1};{1};{1};{1};{1};{1};{1};{1}\n".\
        format(nom,'-')

    #End of the function
    print(best)
    return(best)

def best_blastn_results_all(liste,outputdir,fasta_dir,fasta_extension,database):
    #Container creation : list which will contain all the strains name of the project
    travail = []
    #Reading the name of the strains listed into a .txt file
    with open(liste, "r") as filin:
        for nom in filin:
            travail.append(nom[:-1])

    #Information
    print("\n")
    print("Voici les souches sur lesquelles vous allez travailler {}\n".format(travail))
    print("Le nombre de souches est de {}\n".format(len(travail)))
    print("####################################################################")
    print("############## launch the blastn best gene script ##################")
    print("####################################################################")

    #For each strain in the travai list, launch blastn & blast_nt_best_result functions
    with open("{}{}.csv".format(outputdir,filename),"w") as filout:
        filout.write("Souche;Gene;Position in the genome; Length of the subject;\
        Length of the blasted sequence;ID percentage;Coverage percentage;hsp score;\
        hsp bits;Remarques; Substitutions nb;Nucleotidiques substitutions;Sequence\n")
        for nom in travail:
            #Information
            print("---> {}{}".format(nom,fasta_extension))
            blastn(nom,outputdir,fasta_dir,fasta_extension,database)
            filout.write(blast_nt_best_result(nom,outputdir, threshold))

def clean_xml_blast(path):
    #Listing all the files which corresponding to the research which are in the directory
    for fichier in glob.glob("{}*_blast_xml.txt".format(path)):
        # Cleaning all the xml_files
        subprocess.run(["rm", fichier])

if __name__ == "__main__":

    # Parse command line options
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

    liste = '/Users/Francois/Documents/projets/ecloacae/patric_db/patric_fasta.txt'
    fasta_dir = '/Users/Francois/Documents/projets/ecloacae/patric_db/fasta_patric/'
    database = '/Users/Francois/blast_data_base/atb_resistance/fosfomycin.fsa'
    outputdir = '/Users/Francois/Documents/projets/ecloacae/patric_db/analyses/blast/fosfo/'
    fasta_extension = '.fna'
    threshold = '90'
    filename = 'patric_fasta_fosfo'

    best_blastn_results_all(liste,outputdir,fasta_dir,fasta_extension,database)
    clean_xml_blast(outputdir)
