#!/usr/bin/env python3
## Aout 2018
### Gravey François

### module Loading
import subprocess
import re
import os

mcr_genes = []
regex = re.compile("^>")
with open("/Users/Francois/cge_data_bases/resfinder_db/colistin.fsa", 'r') as filin:
    for ligne in filin:
        if regex.search(ligne):
            print(ligne)
            mcr_genes.append(ligne.split("_")[2][:-1])

print(mcr_genes)
print(len(mcr_genes))
