#!/usr/bin/env python3
## Aout 2018
### Gravey François

### module Loading
import subprocess
import re
import os

def absoluteFilePaths(directory):
   for dirpath,_,filenames in os.walk(directory):
       for f in filenames:
           yield os.path.abspath(os.path.join(dirpath, f))

absoluteFilePaths("/Users/Francois/cge_data_bases/resfinder_db")
for f in os.listdir("/Users/Francois/cge_data_bases/resfinder_db"):
    if f.endswith('.fsa'):
        print(f)
        regex = re.compile("[ATCG]+S[ATCG]+")
        with open("/Users/Francois/cge_data_bases/resfinder_db/{}".format(f), "r") as filin:
            for l in filin:
                if regex.search(l):
                    print(l)
