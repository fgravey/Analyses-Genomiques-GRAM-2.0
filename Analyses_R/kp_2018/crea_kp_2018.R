## Function definitions
plasmides_asso = function(x){
  return(paste(sort(plasmide$Plasmide[grepl(x, plasmide$Strain)]), collapse = "_"))
}

#### Reading files:
epi = read.table("/Volumes/Maxtor/Back_up_pasteur/kp_caen_2018/Data/Kpn_BLSE_rea_2018.csv", 
                 header = TRUE, sep = ";", stringsAsFactors = FALSE)

mlst = read.table("/Volumes/Maxtor/Back_up_pasteur/kp_caen_2018/analyses/MLST/mlst.csv", 
                  header = TRUE, sep = ";", stringsAsFactors = FALSE)

plasmide = read.table("/Volumes/Maxtor/Back_up_pasteur/kp_caen_2018/analyses/plasmidfinder/plasmidfinder_results.csv",
                      header = TRUE, sep = ";", stringsAsFactors = FALSE)

atb = read.csv("/Volumes/Maxtor/Back_up_pasteur/kp_caen_2018/analyses/resistome/Kp_Caen_2018_rea.csv",
                 header = TRUE, sep = ";", stringsAsFactors = FALSE)
pheno = read.table("/Volumes/Maxtor/Back_up_pasteur/kp_caen_2018/Data/Kpn_BLSE_rea_2018.csv", 
                   header = TRUE, sep = ";", stringsAsFactors = FALSE)
pheno = pheno[,-(1:15)]

##### Mooving to the working directory
setwd("/Volumes/Maxtor/Back_up_pasteur/kp_caen_2018/analyses")

## keeping only relevant informations
epi = epi[,-c(7,12,13,15:37)]
epi$N = gsub("KPN-","KpKP", epi$N)
atb = atb[,-c(5,8,9,12:19,21:23,25,27)]

##Merging informations
kp = merge(epi,mlst, by.x = "N", by.y = "Strain")
kp=merge(kp, atb, by.x = "N", by.y = "Souche")

### Creation of a plasmid column into the kp df
plasm = c(paste(sort(plasmide$Plasmide[grepl("KpKP1", plasmide$Strain)]), collapse = "_"),
paste(sort(plasmide$Plasmide[grepl("KpKP2", plasmide$Strain)]), collapse = "_"),
paste(sort(plasmide$Plasmide[grepl("KpKP3", plasmide$Strain)]), collapse = "_"),
paste(sort(plasmide$Plasmide[grepl("KpKP4", plasmide$Strain)]), collapse = "_"),
paste(sort(plasmide$Plasmide[grepl("KpKP5", plasmide$Strain)]), collapse = "_"),
paste(sort(plasmide$Plasmide[grepl("KpKP6", plasmide$Strain)]), collapse = "_"),
paste(sort(plasmide$Plasmide[grepl("KpKP7", plasmide$Strain)]), collapse = "_"),
paste(sort(plasmide$Plasmide[grepl("KpKP8", plasmide$Strain)]), collapse = "_"))
kp$plasmide = plasm

plasm = as.vector(sapply(unique(plasmide$Strain), plasmides_asso))
kp$plasmides = plasm

## wrinting the file
#write.table(kp, "kp_2018_chu_caen_resume.csv", sep = ";", row.names = FALSE)


