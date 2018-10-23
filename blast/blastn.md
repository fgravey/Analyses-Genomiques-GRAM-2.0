# Présentation du pipeline `blastn.py` et `ntoprot.py`

## Gravey François GRAM 2.0

# Présentation du script **Script `blastn.py`**

## Pourquoi ce script ?!

Ce script a été developpé afin de rechercher des gènes d'intérets dans les génomes
bactériens.

## En quel langage est codé ce script ?
Ce script a été developpé en utilisant Python v3.6.0

## Est-ce que ce script a des dépendances ?
**Oui**

  1° - Nécessité d'avoir sur son ordinateur l'excecutable [blast](https://www.ncbi.nlm.nih.gov/books/NBK279690/) qui **doit être ajouté au Path de la machine**.

  2° - Nécéssité d'avoir installé la bibliothèque [Biopython](https://biopython.org)

## Préparation de la base de données avant le makeblastdb :

Afin de pouvoir réaliser le blastn, les gènes d'intérêts doivent être consignés dans un fichier `.fasta`.
De plus **le header du fichier fasta aura une importance considérable***.
Le header devra être écrit de la façon ci-contre `>Nomdugene_informationscomplementaires-informationsnongardées. En effet, le script va rechercher le nom du gène pour écrire les séquences nucléotidiques, ainsi, **TOUS** les header doivent commencer par cette information. De plus, les informations contenues dans le header qui seront après le back slash ne seront pas gardées pour la suite du travail. Par conséquent, faites **très ATTENTION** à la façon de préparer votre fichier.

Exemples :

Pour un header définit comme suit : `>ampC_ECL_34978` le gène qui sera extrait par le script sera le suivant : **ampC**.

Pour un header définit comme suit : `>mcmI_Microcin_M_Strain_CA46-accession : AJ515251.1`, le gène considéré par le script sera le suivant : **mcmI_Microcin_M_Strain_CA46**, aussi, c'est à l'utilisateur de définir les informations qu'il souhaite associer, ou non, au gène sur lequel il travail.

Le deuxième exemple peut être utile, lorsque vous travaillez à partir de séquences proches pour un même gène. Dans ce cas, il est intéressant, de "garder", la souche dont la séquence a été récupérée.

## Quels sont les inputs du script ?
  Pour fonctionner le script réclame sept inputs :

  1°- `-l` est un fichier au format **.txt** qui contient la liste de l'ensemble des préfixes des fichiers fasta sur lesquels vous souhaitez travailler. Vous devez renseigner le **chemin ABSOLU** du fichier ;

  2°- `-f` est le **chemin ABSOLU** du repertoire qui contient les fichiers fasta dans lesquels vous souhaitez rechercher les gènes d'intérêt ;

  3°- `-o` le **chemin ABSOLU** du repertoire dans lequel vous souhaitez que les résultats des blasts prennent place  ;

  4°- `-db` le **chemin ABSOLU** vers le fichier qui a servi de référence à la base de données blast **ATTENTION** ce fichier soit être indexé selon la commande `makeblastdb`, plus d'informations sont disponibles sur ce lien [makeblastdb](https://www.ncbi.nlm.nih.gov/books/NBK279688/) ;

  5°- `-e` est **l'extension** des fichiers fasta par exemple `.fasta` ou `.agp.fasta` ou `.awked.fasta` ;

  6°- `filename`est le nom que vous souhaitez donner au fichier de résultat par défaut, le nom du fichier sera **summary_blast** ;

  7°- `-t` est le seuil de coverage que vous souhaitez appliquer aux résultats des blasts, par défaut, la valeur est **80**.


## Quels sont les outputs de ce scripts ?
  Les outputs du scripts sont **NOMBREUX** :

  1°- Pour chaque souche considérée, un répertoire de travail portant le nom de la souche va automatiquement se créer : `outputdir/**nomdelasouche**`. Ce dossier va contenir différents fichiers :

  - Un fichier `nomdelasouche_blast_xml.txt` qui sera parsé par [Biopython](https://biopython.org) afin d'en extraire les informations utiles,

  - Un fichier `nomdelasouche_blast.txt` qui permettra à l'utilisateur de pouvoir regarder dans un format facile, les résultats des blasts,

  - Un **dossier** `outputdir/**nomdelasouche**/nomdelasouche_genes_sequences_fasta` qui contient les séquences nucléotidiques correspondant aux gènes d'intérêts extraites du génome. Pour chaque gène retrouvé, un fichier fasta dédié sera créé et nommé sous la forme suivante : `gene_nomdelasouche_nt_sequence.fasta`. Le header du fichier `.fasta` sera alors le suivant `>nomdelasouche_gene_sequence`. Le header peut contenir l'information **reverse** qui indique que le brin du génome assemblé, n'est pas le brin codant du gène recherché,

  2°- Un dossier sera automatiquement créé, nommé `outputdir/**multifasta_nt**`. Ce dernier renfermera différents fichiers `.fasta`; un par gène d'intérêt retrouvé au moins une fois. Chaque fichier `.fasta`, renferme la séquence de référence du gène d'intérêt, ainsi que l'ensemble des gènes retrouvés dans les différents génomes testés. Par conséquent, le script créé un fichier multifasta qui est prêt à l'emploi pour réaliser un alignment de séquences via l'outil [Muscle](https://www.ebi.ac.uk/Tools/msa/muscle/) par exemple.

  3°- Un fichier au format `.csv` nommé `nomdefichierchoisi_per_strain_blastn.csv` qui contient différentes informations à propos des gènes retrouvés pour chaque génome inclus dans le blast. Les informations disponibles seront les suivantes :

  - Nom de la souche,
  - Nom du gene,
  - Contig dans lequel le gène a été retrouvé,
  - Longueur du gène de référence,
  - Longueur de la séquence nucléotidique dans le génome correspondant au gène recherché,
  - Pourcentage de coverage,
  - Pourcentage d'identité,
  - Position du gène dans le contig,
  - Sens du brain d'ADN sur lequel le gène a été retrouvé
  - Remaque concernant de possibles pertes en 3' ou en 5' de séquence,
  - Nombre de substitutions nucleotidiques,
  - Position des substitutions nucleotidiques.

  4°- Un fichier au format `.csv` nommé `nomdefichierchoisi_per_gene_blastn.csv` qui contient différentes informations à propos des gènes retrouvés pour chaque génome inclus dans le blast. Les informations disponibles seront les suivantes :

  - Nom de la souche,
  - Nom du gene,
  - Contig dans lequel le gène a été retrouvé,
  - Longueur du gène de référence,
  - Longueur de la séquence nucléotidique dans le génome correspondant au gène recherché,
  - Pourcentage de coverage,
  - Pourcentage d'identité,
  - Position du gène dans le contig,
  - Sens du brain d'ADN sur lequel le gène a été retrouvé
  - Remaque concernant de possibles pertes en 3' ou en 5' de séquence,
  - Nombre de substitutions nucleotidiques,
  - Position des substitutions nucleotidiques.

  Les différences entre les fichiers `nomdefichierchoisi_per_strain_blastn.csv et `nomdefichierchoisi_per_gene_blastn.csv` ne sont que de pure forme, les deux fichiers contiennent les même informations. Le premier centré sur chaque souche, le second centré sur chaque gène.

# Présentation du script **Script `ntoprot.py`**

## Pourquoi ce script ?!

Ce script a été developpé afin de rechercher des gènes d'intérets dans les génomes
bactériens.

## En quel langage est codé ce script ?
Ce script a été developpé en utilisant Python v3.6.0

## Est-ce que ce script a des dépendances ?
**Oui**

  1° - Nécessité d'avoir sur son ordinateur l'excecutable [blast](https://www.ncbi.nlm.nih.gov/books/NBK279690/) qui **doit être ajouté au Path de la machine**.

  2° - Nécéssité d'avoir installé la bibliothèque [Biopython](https://biopython.org)
