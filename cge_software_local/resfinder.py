#!/usr/bin/env python3
## Aout 2018
### Gravey François

### module Loading
import subprocess
import re

## function definition
def indices(mylist, value):
    """ Return the indices of a specific value contained in a list"""
    return [i for i,x in enumerate(mylist) if x==value]

def format_atb(atb):
    """Change special character in order to look for it using regex and
    the re. module"""
    atb = re.sub(r"\(",r"\\(",atb)
    atb = re.sub(r"\)",r"\\)",atb)
    atb = re.sub(r"\[",r"\\[",atb)
    atb = re.sub(r"\]",r"\\]",atb)
    return atb

def looking_for_missmatch(gene,nom):
    """ Look for nucleotide substitution between the reference genome and the
    genome of the strain. Using the blast results presented in the results.txt
    file.
        - Inputs : name of the gene & name of the strain
        - Output : list which contains alerte if there are some nucleotide
        substitutions"""

    ## Path definition
    outputdir = "{}/{}".format(outdir,nom)

    ## definition of the four regular expressions
    regex = re.compile("{}:".format(format_atb(gene))) #looking for the gene in the results files
    regex_1 = re.compile("Resistance gene seq: ") #looking for the sequence of the reference genome
    regex_2 = re.compile("Hit in genome:       ") #lookng for the sequence of the strain
    regex_3 = re.compile(" "*21) #looking for the adaquation line between the two genomes
    #regex_4 = re.compile('--------------------------------------------------------------------------------')
    regex_4 = re.compile('-'*80) #looking for the end of the result

    ## definition of three containers : souche which will contain the
    # nucleotide sequence of the strain, genome_ref which will contain
    # the nucleotide sequence of the reference genome, and match which
    # will contain the all the match (perfect or not) of the two
    # sequences
    souche = []
    genome_ref = []
    match = []

    with open("{}/results.txt".format(outputdir), 'r') as filin:
        line = filin.readline()
        while line:
            if regex.search(line): #parsing the file until find the gene of interest
                line = filin.readline() # then for all the following lines we look for
                # the reference genome or the genome of the strain or the match between
                # the two genomes
                while regex_4.search(line) == None:
                    if regex_1.search(line):
                        seq = regex_1.sub('', line)
                        genome_ref.append(seq[:-1])
                    elif regex_2.search(line):
                        target = regex_2.sub('', line)
                        souche.append(target[:-1])
                    elif regex_3.search(line):
                        space = regex_3.sub('', line)
                        match.append(space[:-1])
                    line = filin.readline()
                    #if regex_4.search(line) != None:
                        #print("Terminé")
                break #When the regex_4 is != None the loops stop

            else:
                line=filin.readline()

    ## remooving coma form the containers
    genome_ref = "".join(genome_ref)
    souche = "".join(souche)
    match = "".join(match)

    ## looking for unperfect match in the macth container illustrated
    ## by a blanck space
    alerte = []
    sub = indices(match, " ")
    if sub:
        for sb in sub:
            if not alerte:
                alerte.append("ATTENTION, il y a {} substitution(s): {} --> {} en position {}".format\
            (len(sub),genome_ref[sb], souche[sb], sb))
            else:
                alerte.append(" {} --> {} en position {}".format\
            (genome_ref[sb], souche[sb], sb))
        alerte = ",".join(alerte)
    else:
        alerte = "-"

    #print(genome_ref)
    #print(match)
    #print(souche)
    #print(alerte)
    return(alerte)

def results_resfinder(outdir,nom):

    #### unchanging path
    resfinder = '/Users/Francois/cge_softwares/resfinder/resfinder.pl'
    database = '/Users/Francois/cge_data_bases/resfinder_db'
    blast_path = "/Users/Francois/cge_softwares/blast-2.2.26"

    #path definition and directories creation
    outputdir = "{}/{}".format(outdir,nom)
    outputdir_blast = "{}/blast".format(outputdir)
    subprocess.run(["mkdir", outputdir])
    subprocess.run(["mkdir", outputdir_blast])

    ## definition of containers
    header = []
    antibio = []

    molecules = ["beta-lactam","aminoglycoside","quinolone", "colistin", \
     "fosfomycin","fusidicacid", "macrolide", "nitroimidazole", "oxazolidinone",\
    "phenicol", "rifampicin", "sulphonamide", "tetracycline","trimethoprim"]
    #molecules = ["aminoglycoside","macrolide","colistin", "fosfomycin",]
    for atb in molecules:

    # for each fasta in the list "travail" and for each antibiotic in the list
    # we launch the resfinder program
        fasta = "{}/{}{}".format(fasta_dir, nom, fasta_extension)
        print("##################################################")
        print("working on {}{} file {}".format(nom, fasta_extension, atb))
        print("##################################################")
        subprocess.run(["perl", "{}".format(resfinder), "-d", "{}".format(database), \
        "-i", "{}".format(fasta), "-o", "{}".format(outputdir), "-k", "95.00", "-l",\
         "0.60", "-b", "{}".format(blast_path), "-a", "{}".format(str(atb))])

        ## parsing the results_table.txt file and adding all the lines into a
        ## container called resultat
        with open("{}/results_table.txt".format(outputdir), 'r') as filin:
            resultat = []
            for l in filin:
                resultat.append(l[:-1])

            with open("{}/resfinder_results_{}_{}.txt".format(outputdir,atb,nom), 'w') as filout:
                filout.write("{}\t{}\t{}\n".format("Souche",resultat[1],"Remarques"))
                for i in range(2,len(resultat)-1):
                    filout.write("{}\t{}\t{}\n".format(nom,resultat[i],\
                    looking_for_missmatch(resultat[i].split("\t")[0],nom)))

        if len(resultat) == 3:
            if molecules.index(atb) == 0:
                header.append("{}\t{}\t{}\t".format("Souche",atb,"Remarques"))
                antibio.append("{}\t{}\t{}\t".format(nom,"-","-"))

            elif int(molecules.index(atb)) == int((len(molecules)-1)):
                header.append("{}\t{}\n".format(atb,"Remarques"))
                antibio.append("{}\t{}\n".format("-","-"))
            else:
                header.append("{}\t{}\t".format(atb,"Remarques"))
                antibio.append("{}\t{}\t".format("-","-"))
            print(header)
            print(antibio)
        if len(resultat) > 3:
            if molecules.index(atb) == 0:
                header.append("{}\t{}\t{}\t".format("Souche",atb,"Remarques"))
                genes = []
                remarque = []
                for i in range(2,len(resultat)-1):
                    genes.append("{}_({}%:{})".format(resultat[i].split("\t")[0], resultat[i].split("\t")[1],\
                    resultat[i].split("\t")[2]))
                    remarque.append("{}".format(looking_for_missmatch(resultat[i].split("\t")[0],nom)))
                genes = "_".join(genes)
                remarque = " | ".join(remarque)
                antibio.append("{}\t{}\t{}\t".format(nom,genes,remarque))

            elif int(molecules.index(atb)) == int((len(molecules)-1)):
                header.append("{}\t{}\n".format(atb,"Remarques"))
                genes = []
                remarque = []
                for i in range(2,len(resultat)-1):
                    genes.append("{}_({}%:{})".format(resultat[i].split("\t")[0], resultat[i].split("\t")[1],\
                    resultat[i].split("\t")[2]))
                    remarque.append("{}".format(looking_for_missmatch(resultat[i].split("\t")[0],nom)))
                genes = "_".join(genes)
                remarque = " | ".join(remarque)
                antibio.append("{}\t{}\n".format(genes,remarque))
            else:
                header.append("{}\t{}\t".format(atb,"Remarques"))
                genes = []
                remarque = []
                for i in range(2,len(resultat)-1):
                    genes.append("{}_({}%:{})".format(resultat[i].split("\t")[0], resultat[i].split("\t")[1],\
                    resultat[i].split("\t")[2]))
                    remarque.append("{}".format(looking_for_missmatch(resultat[i].split("\t")[0],nom)))
                genes = "_".join(genes)
                remarque = " | ".join(remarque)
                antibio.append("{}\t{}\t".format(genes,remarque))

        ## copying the results.txt file and renaming it in order to keep a copy
        ## of it for each antibiotic and each strain
        with open("{}/results.txt".format(outputdir), 'r') as filin:
            with open("{}/results_blast_{}_{}.txt".format(outputdir_blast, atb, nom), 'w') as filout:
                for l in filin:
                    filout.write(l)

    header = "".join(header)
    antibio = "".join(antibio)
    final = []
    final.append(header)
    final.append(antibio)
    return(final)


##### changing path
liste = '/Volumes/Maxtor/lugdunensis/liste_souches_lugdu.txt'
#working list which contains all the name of the working files
fasta_dir = '/Volumes/Maxtor/lugdunensis/fasta_assembled'
#directory which contains the fasta
outdir = '/Volumes/Maxtor/lugdunensis/analyses/resistome' #path to output files
fasta_extension = '.scfd.fasta'

#listing all the files which will be working on
travail = []
with open(liste, 'r') as filin:
    for nom in filin:
        travail.append(nom[:-1])

print("Voici la liste des fichiers sur lesquels vous allez travailler", travail)
print("Le nombre de fichier est de : {}".format(len(travail)))

#On parcours la liste et pour chaque nom, on execute resfinder pour tous les
#atb de la liste
for nom in travail:
    outputdir = "{}/{}".format(outdir,nom)
    if travail.index(nom) == 0:
        with open("{}/summary_resfinder.txt".format(outdir), 'w') as filout:
            for l in results_resfinder(outdir,nom):
                filout.write(l)
    else:
        resultat = []
        with open("{}/summary_resfinder.txt".format(outdir), 'a') as filout:
            for l in results_resfinder(outdir,nom):
                resultat.append(l)
            filout.write(resultat[1])
    print("##################################################")
    print("Cleaning results tempory files ")
    print("##################################################")
    subprocess.run(["rm", "-rf", "{}/results.txt".format(outputdir)])
    subprocess.run(["rm", "-rf", "{}/Hit_in_genome_seq.fsa".format(outputdir)])
    subprocess.run(["rm", "-rf", "{}/Resistance_gene_seq.fsa".format(outputdir)])
    subprocess.run(["rm", "-rf", "{}/results_tab.txt".format(outputdir)])
    subprocess.run(["rm", "-rf", "{}/results_table.txt".format(outputdir)])



#print(results_resfinder(outdir,'C15-10'))
