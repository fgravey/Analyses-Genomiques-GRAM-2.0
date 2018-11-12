#!/usr/bin/env python3
## Octobre 2018
### Gravey François

### module Loading
import subprocess
import re
import json

##### changing path
liste = '/Volumes/Maxtor/lugdunensis/liste_souches_lugdu_test.txt'
#working list which contains all the name of the working files
fasta_dir = '/Users/Francois/Documents/projets/lugdu/fasta_assembled/'
#directory which contains the fasta
outputdir = '/Users/Francois/Documents/projets/lugdu/analyses/virulencefinder/' #path to output files
fasta_extension = '.scfd.fasta'
nom = '16491'
specie = 'staphylococcus'
filename = 'essai'



def virulencefinder_strain(nom,fasta_dir,fasta_extension,outputdir,specie):
    #### unchanging path
    virulencefinder = '/Users/Francois/Desktop/virulencefinder/virulencefinder.py'
    database = '/Users/Francois/Desktop/virulencefinder_db/'

    #Variable definition
    fasta = "{}{}{}".format(fasta_dir,nom,fasta_extension)
    res = [] #list which will contains all the virulence genes recovered

    #Definition databases:
    if specie == 'staphylococcus':
        db = ['s.aureus_exoenzyme', 's.aureus_hostimm', 's.aureus_toxin']

    elif specie == 'ecoli':
        db = ['virulence_ecoli','stx']

    elif specie == 'enterococcus':
        db = ['virulence_ent']

    elif specie == 'listeria':
        db = ['listeria']
    print(db)
    for research in db:
        #Launching the virulencefinder.py from cge script
        subprocess.run(["python", "{}".format(virulencefinder), "-i", "{}".format(fasta),\
        "-o","{}".format(outputdir),"-tmp","{}".format(outputdir), "-mp", "blastn",\
        "-p", "{}".format(database), "-x", "-d", "{}".format(research)])

        #Parsing the virulencefinder results in the .json file and add all the information in the data variable
        with open('{}data.json'.format(outputdir), 'r') as filin:
            data = json.load(filin)

        #Varibale definition:
        donnees = data['virulencefinder']['results']

        for espece in donnees.keys():
            for base_donnee in donnees[espece].keys():
                if donnees[espece][base_donnee] == 'No hit found':
                    res.append("{0};{1};No hit found;{2};{2};{2};{2};{2}\n"\
                    .format(nom,base_donnee,'-'))
                else:
                    for gene in donnees[espece][base_donnee].keys():
                        resultat = donnees[espece][base_donnee][gene]
                        res.append("{};{};{};{};{};{};{};{}\n".format(nom,\
                        base_donnee,resultat['virulence_gene'],\
                        resultat['identity'],resultat['coverage'],\
                        resultat['contig_name'],resultat['positions_in_contig']\
                        ,resultat['protein_function']))
    #End of the function
    return(res)

def virulencefinder_strain(nom,fasta_dir,fasta_extension,outputdir,specie):
    listing all the files which will be working on
    travail = []
    with open(liste, 'r') as filin:
        for nom in filin:
            travail.append(nom[:-1])

    # Information
    print("\n")
    print("Voici la liste des fichiers sur lesquels vous allez travailler : {}\n"\
    .format("".join(travail)))
    print("Le nombre de fichier est de : {}\n".format(len(travail)))
    print("#########################################################")
    print("############ Launching serotypefinder ###################")
    print("#########################################################")

    with open("{}{}".format(outputdir,filename),'w') as filout:
        # Information
        print("---> {}{}.csv".format(nom,extension))

        # Writing the output file
        filout.write("Souche;database;Gene;Identity;Coverage;Contig;Position in contig;Protein function\n")
        #for each strain in the list, virulencefinder_strain is executed
        for nom in travail:
            filout.write(virulencefinder_strain(nom,fasta_dir,fasta_extension,outputdir,specie))
