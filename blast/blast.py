#!/usr/bin/env python3
## Aout 2018
### Gravey François

### module Loading
import subprocess
import re
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML

##Function defintion
def blastn(nom,outputdir,fasta_dir):
    """Make a directory specific to each strain wich contains two blast files
    first one is a .txt file and the second one is an .xml file.
    The xml file will be parse by bipython
    Inputs : nom = name of the strain
             outputdir = directory which will contains all the subdirectories
             fasta_dir = directory which contains the fasta files

    Outputs : a specific output by strain which contains the two fasta files
    """
    outputdir = "{}/{}".format(outputdir,nom)
    subprocess.run(["mkdir", outputdir])
    #launch blast between all the files in the list and the database
    fasta = "{}/{}.{}".format(fasta_dir,nom,fasta_extension)
    subprocess.run(["blastn", "-query", "{}".format(fasta), "-db", \
    "{}".format(database), "-out", "{}/{}_blast_xml.txt".format(outputdir,nom),\
    '-outfmt', '5'])
    subprocess.run(["blastn", "-query", "{}".format(fasta), "-db", \
    "{}".format(database), "-out", "{}/{}_blast.txt".format(outputdir,nom)])


def blast_nt_result(nom,outputdir):
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
    outputdir = "{}/{}".format(outputdir,nom)
    with open("{}/{}_blast_xml.txt".format(outputdir,nom)) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        E_VALUE_THRESH = 0.01
        compteur = 0 #count the number of genes found in the strain's genome
        resultats = []
        for blast_record in blast_records:
            #print(dir(blast_record))
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
                        compteur = compteur + 1
                        resultats.append("{};{};{};{};{};{};{};{}\n".format(souche,\
                        gene, contig,longueur, id, positive, gaps, diff))
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

##### changing path
liste = '/Volumes/Maxtor/Back_up_pasteur/E_cloacae_F_Guerrin/ecloacae.txt'
#working list which contains all the name of the working files
fasta_dir = '/Volumes/Maxtor/Back_up_pasteur/E_cloacae_F_Guerrin/raw/fasta/agp_fasta'
#directory which contains the fasta
outputdir = '/Volumes/Maxtor/Back_up_pasteur/E_cloacae_F_Guerrin/analyses/blast' #path to output files
fasta_extension = 'agp.fasta'
database = '/Users/Francois/blast_data_base/ecloacae/interet/genes_ecc_fg.fasta'
nom_fichier = 'ecc_fg_FAY_blast'

with open("{}/{}.csv".format(outputdir,nom_fichier), 'w') as filout:
    for results in blast_nt_result_filout(liste, fasta_dir, outputdir, fasta_extension, database):
        filout.write(results)
