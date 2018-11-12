# Notes à propos de la gestion des bases de données cge pour fonctionnement des scripts

## Gravey François GRAM 2.0
### Octobre 2018

## mlst_db, virulencefinder_db et serotypefinder_db
Les bases de données mlst_db_cge sont conçues pour fonctionner à la fois avec `kma` mais aussi avec `blastn`.
Pour indexer les bases de données `kma`, utiliser la commande suivante : `python3 INSTALL.py kma_index`

Pour indexer la base de données `blastn`, de mlst_db, un script python nommé index_blastn_db.py a été spécialement créé.

Pour indexer les bases de données **virulencefinder_db** et **serotypefinder_db**, le plus simple est de lancer la commande bash suivante : `for f in *.fsa; do makeblastdb -in $f -dbtype 'nucl'; done`. Ainsi, tous les fichiers `.fsa`contenus dans le répertoire seront indexés et prêts à être utilisés par `blastn`.


## Liste des bases de données cge utilisés:

- fimtyper --> script python ok
- mlst
- plasmidfinder --> script python ok
- pmlst --> script python ok
- resfinder --> script python ok
- serotypefinder updated 2018-08-03
- virulencefinder version 2018-10-12
