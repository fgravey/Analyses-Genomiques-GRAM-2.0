## Loading functions
source("/Users/Francois/Documents/These_de_science/scripts/Analyses_R/function.R")

#### Reading informations files EPI
## rea project
epi_rea = read.csv("/Volumes/Maxtor/Back_up_pasteur/Rea_neonat/Data/Clinique/clinique_rea.csv",
                   header = TRUE, sep = ";", stringsAsFactors = FALSE)

epi_rea = epi_rea[,c("Nom.etude","Patient", "Ne.e..le", "Sexe", "Date.entree", "Date.de.prel.", "Nature")]

### 2018 E. coli ESBL project
epi_2018 = read.csv("/Volumes/Maxtor/Back_up_pasteur/blse_ecoli_2018_ajout/data/epi/listing_souches_nouvelles.csv",
                    header = TRUE, sep = ";", stringsAsFactors = FALSE)

epi_2018 = epi_2018[,c("N", "Patient", "Ne.e..le", "Sexe", "Date.entree", "Date.de.prel.", "Nature")]
epi_2018$Nom.etude = gsub("ECO-", "", epi_2018$Nom.etude)
epi_2018$Nom.etude = gsub("-", "", epi_2018$Nom.etude)

### 2012 2017 E. coli ESBL project
epi_1217 = read.csv("/Volumes/Maxtor/Back_up_pasteur/BLSE_Pasteur/data/donnees_cliniques/epi.csv",
                    header = TRUE, sep = ";", stringsAsFactors = FALSE)

epi_1217 = epi_1217[,c("Nom.etude", "Patient", "Ne.e..le", "Sexe", "Date.entree", "Date.de.prel.", "Nature")]
epi_1217$Nom.etude = gsub("B", "C", epi_1217$Nom.etude)

## Uniformization of the colnames
colnames(epi_2018) = c("Nom.etude", "Patient", "Ne.e..le", "Sexe", "Date.entree", "Date.de.prel.", "Nature")

#### Reading informations files Resistome using atb function
atb_1217 = atb("/Volumes/Maxtor/Back_up_pasteur/BLSE_Pasteur/Analyses/Resistance/blse_1217_Resfinder_260218.txt")
atb_rea = atb("/Volumes/Maxtor/Back_up_pasteur/Rea_neonat/Analyses/Resistance/rea_Resfinder.txt")
atb_2018 = atb("/Volumes/Maxtor/Back_up_pasteur/blse_ecoli_2018_ajout/analyses/Resistance/2018_Resfinder.txt")

#### Reading informations files MLST 
mlst_2018 = mlst("/Volumes/Maxtor/Back_up_pasteur/blse_ecoli_2018_ajout/analyses/MLST/2018_MLSTfinder.txt")
mlst_rea = mlst("/Volumes/shigella-ngs/EcCaen/Rea_neonat/Analyses/MLST/rea_MLSTfinder.txt")
