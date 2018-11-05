#!/usr/bin/env python3
## Novembre 2018
### Gravey François

import subprocess
from blast import blastn
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML


nom = 'FAY1'
fasta_dir = '/Volumes/Maxtor/Back_up_pasteur/E_cloacae_F_Guerrin/raw/fasta/agp_fasta/'
database = '/Users/Francois/blast_data_base/ecloacae/clusters/clusters_ecloacae.fsa'
outputdir = '/Users/Francois/Desktop/essai/'
fasta_extension = '.agp.fasta'
threshold = '90'


#blastn(nom,outputdir,fasta_dir,fasta_extension,database)

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
    #multifasta_nt_sequence_dir = "{}/{}_genes_sequences_fasta".format(inputdir,nom)
    #subprocess.run(["mkdir", multifasta_nt_sequence_dir])

    ## Reading the xml format of the blast results file
    with open("{}/{}_blast_xml.txt".format(inputdir,nom)) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        E_VALUE_THRESH = 0.0000001
        resultats = []
        best = []
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

                            ###looking for 3' or 5' deletions
                            if sbjct_start != 1 and sbjct_end != 1:
                                remarque.append("tronquee en 5'")
                            if sbjct_start != subject_longueur and sbjct_end != subject_longueur:
                                remarque.append("tronquee en 3'")
                            if not remarque :
                                remarque.append("-")

                            nb_substitutions = len(substitutions)
                            if not substitutions:
                                substitutions.append('-')

                            if not best:
                                best = "{};{};{};{};{};{}".format(nom,gene,contig,\
                                subject_longueur, alignement_longueur,\
                                perc_ident,perc_coverage,hsp.score,hsp.bits,\
                                nb_substitutions)

                            if float(best[2]) > float(hsp.bits):
                                print("il y a mieux")
                            else:
                                print(gene)


print(blast_nt_result(nom,outputdir, threshold))
