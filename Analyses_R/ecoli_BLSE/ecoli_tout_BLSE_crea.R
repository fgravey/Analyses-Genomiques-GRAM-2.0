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
## Uniformization of the colnames
colnames(epi_2018) = c("Nom.etude", "Patient", "Ne.e..le", "Sexe", "Date.entree", "Date.de.prel.", "Nature")
epi_2018$Nom.etude = gsub("ECO-", "", epi_2018$Nom.etude)
epi_2018$Nom.etude = gsub("-", "", epi_2018$Nom.etude)

### 2012 2017 E. coli ESBL project
epi_1217 = read.csv("/Volumes/Maxtor/Back_up_pasteur/BLSE_Pasteur/data/donnees_cliniques/epi.csv",
                    header = TRUE, sep = ";", stringsAsFactors = FALSE)

epi_1217 = epi_1217[,c("Nom.etude", "Patient", "Ne.e..le", "Sexe", "Date.entree", "Date.de.prel.", "Nature")]
epi_1217$Nom.etude = gsub("B", "C", epi_1217$Nom.etude)



#### Reading informations files Resistome using atb function
atb_1217 = atb("/Volumes/Maxtor/Back_up_pasteur/BLSE_Pasteur/Analyses/Resistance/blse_1217_Resfinder_260218.txt")
atb_rea = atb("/Volumes/Maxtor/Back_up_pasteur/Rea_neonat/Analyses/Resistance/rea_Resfinder.txt")
atb_2018 = atb("/Volumes/Maxtor/Back_up_pasteur/blse_ecoli_2018_ajout/analyses/Resistance/2018_Resfinder.txt")

#### Reading informations files MLST 
mlst_2018 = mlst("/Volumes/Maxtor/Back_up_pasteur/blse_ecoli_2018_ajout/analyses/MLST/2018_MLSTfinder.txt")
mlst_rea = mlst("/Volumes/Maxtor/Back_up_pasteur/Rea_neonat/Analyses/MLST/rea_MLSTfinder.txt")
mlst_1217 = mlst("/Volumes/Maxtor/Back_up_pasteur/BLSE_Pasteur/Analyses/MLST/1217_MLSTfinder.txt")

### Readling the serotype files :
sero_2018 = sero("/Volumes/Maxtor/Back_up_pasteur/blse_ecoli_2018_ajout/analyses/Serotype/2018_SerotypeFinderColi.txt")
sero_rea = sero("/Volumes/Maxtor/Back_up_pasteur/Rea_neonat/Analyses/Serotype/ecoli_rea_SerotypeFinderColi.txt")
sero_1217 = sero("/Volumes/Maxtor/Back_up_pasteur/BLSE_Pasteur/Analyses/Serotype/1217_SerotypeFinderColi.txt")

### Reading the FimH files
fim_rea = fim("/Volumes/Maxtor/Back_up_pasteur/Rea_neonat/Analyses/FimH/rea_fimH_typing.txt")
fim_1217 = fim(path = "/Volumes/Maxtor/Back_up_pasteur/BLSE_Pasteur/Analyses/FimH/1217_resulats_fimH_typing.txt")
fim_2018 = fim("/Volumes/Maxtor/Back_up_pasteur/blse_ecoli_2018_ajout/analyses/FimH/2018_resulats_fimH_typing.txt")

### Reading the Virulence files
viru_1217 = viru("/Volumes/Maxtor/Back_up_pasteur/BLSE_Pasteur/Analyses/Virulence/1217_Virulencefinder.txt")
viru_rea = viru("/Volumes/Maxtor/Back_up_pasteur/Rea_neonat/Analyses/Virulence/rea_virulencefinder.txt")
viru_2018 = viru("/Volumes/Maxtor/Back_up_pasteur/blse_ecoli_2018_ajout/analyses/Virulence/2018_virulencefinder.txt")

coli2018 = merge(epi_2018, mlst_2018, by.x = "Nom.etude", by.y = "file", all.y = TRUE)
coli2018 = merge(coli2018, sero_2018, by.x = "Nom.etude", by.y = "file", all.y = TRUE)
coli2018 = merge(coli2018, fim_2018, by.x = "Nom.etude", by.y = "Souche", all.y = TRUE)
coli2018 = merge(coli2018, viru_2018, by.x = "Nom.etude", by.y = "file", all.y = TRUE)
coli2018 = merge(coli2018, atb_2018, by.x = "Nom.etude", by.y = "File", all.y = TRUE)

colirea = merge(epi_rea, mlst_rea, by.x = "Nom.etude", by.y = "file", all.y = TRUE)
colirea = merge(colirea, sero_rea, by.x = "Nom.etude", by.y = "file", all.y = TRUE)
colirea = merge(colirea, fim_rea, by.x = "Nom.etude", by.y = "Souche", all.y = TRUE)
colirea = merge(colirea, viru_rea, by.x = "Nom.etude", by.y = "file", all.y = TRUE)
colirea = merge(colirea, atb_rea, by.x = "Nom.etude", by.y = "File", all.y = TRUE)

coli1217 = merge(epi_1217, mlst_1217, by.x = "Nom.etude", by.y = "file", all.y = TRUE)
coli1217 = merge(coli1217, sero_1217, by.x = "Nom.etude", by.y = "file", all.y = TRUE)
coli1217 = merge(coli1217, fim_1217, by.x = "Nom.etude", by.y = "Souche", all.y = TRUE)
coli1217 = merge(coli1217, viru_1217, by.x = "Nom.etude", by.y = "file", all.y = TRUE)
coli1217 = merge(coli1217, atb_1217, by.x = "Nom.etude", by.y = "File", all.y = TRUE)

colim2 = rbind(coli1217, colirea)
colim2 = rbind(colim2, coli2018)
colim2 = as.data.frame(t(apply(colim2, 1, clean_percent)))
colim2$ST = gsub("Incomplete_ICD[0-9]+_[0-9]+/[0-9]+", "ND", colim2$ST)
resume_years(colim2$ST, colim2)

setwd("/Volumes/Maxtor/Back_up_pasteur")
write.xlsx(colim2, file = "stats_colim2_avant_2018_essai_2.xlsx",
           sheetName = "tout", row.names = FALSE, append = FALSE)
write.xlsx(resume_years(colim2$ST,colim2), file = "stats_colim2_avant_2018_essai_2.xlsx",
           sheetName = "mlst", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(colim2$OH,colim2), file = "stats_colim2_avant_2018_essai_2.xlsx",
           sheetName = "OH", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(colim2$Fimtype,colim2), file = "stats_colim2_avant_2018_essai_2.xlsx",
           sheetName = "Fimtype", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(colim2$virulence,colim2), file = "stats_colim2_avant_2018_essai_2.xlsx",
           sheetName = "virulence", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(colim2$aminoglycoside,colim2), file = "stats_colim2_avant_2018_essai_2.xlsx",
           sheetName = "aminosides", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(colim2$beta.lactamase,colim2), file = "stats_colim2_avant_2018_essai_2.xlsx",
           sheetName = "beta_lactamases", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(colim2$quinolone,colim2), file = "stats_colim2_avant_2018_essai_2.xlsx",
           sheetName = "quinolones", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(colim2$fosfomycin,colim2), file = "stats_colim2_avant_2018_essai_2.xlsx",
           sheetName = "fosfomycine", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(colim2$trimethoprim,colim2), file = "stats_colim2_avant_2018_essai_2.xlsx",
           sheetName = "trimethoprim", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(colim2$macrolide,colim2), file = "stats_colim2_avant_2018_essai_2.xlsx",
           sheetName = "macrolides", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(colim2$phenicol,colim2), file = "stats_colim2_avant_2018_essai_2.xlsx",
           sheetName = "phenicol", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(colim2$sulphonamide,colim2), file = "stats_colim2_avant_2018_essai_2.xlsx",
           sheetName = "sulphonamide", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(colim2$tetracyclinet,colim2), file = "stats_colim2_avant_2018_essai_2.xlsx",
           sheetName = "tetracycline", row.names = FALSE, append = TRUE)
write.csv()

