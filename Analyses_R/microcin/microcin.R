#### mooving to the working directory
library(xlsx)

setwd("/Volumes/Maxtor/Back_up_pasteur/microcin/analyses")
stec15 = read.table("microcin_Ec2015.txt")
stec16 = read.table("microcin_Ec2016.txt")
blse = read.csv("blse_caen.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE)

### Keeping only few data:
stec15 = stec15[,c(1,2)]
stec16 = stec16[,c(1,2)]
colnames(stec15) = c("souche", "gene")
colnames(stec16) = c("souche", "gene")

blse = blse[,c(1,2)]

final = rbind(blse,stec15)
final = rbind(final,stec16)

final$souche = gsub("_NODE_[0-9]+_covk_[0-9]+\\.[0-9]+_covr_[0-9]+\\.[0-9]+_cutoff_[0-9]+_taxo_[A-Z]+","",final$souche)
final$souche = gsub("_NODE_[0-9]+_covk_[0-9]\\.[0-9]+_covr_[0-9]+\\.[0-9]+_cutoff_[0-9]","",final$souche)
final$souche = gsub("_NODE_[0-9]+_covk_[0-9]+\\.[0-9]+_covr_[0-9]+\\.[0-9]+_cutoff_[0-9]+","",final$souche)
final$souche = gsub(".agp.fasta", "", final$souche)
pc = as.data.frame(sort(table(final$gene), decreasing = TRUE), row.names = NULL)

##nb of strains :
# blse = 56
# rea = 35
# blse ajout = 27
# stec 2015 = 118
# stec 2016 = 186

pc$percentage = round(pc$Freq/(56+35+27+118+186)*100, 2)
colnames(pc) = c("gene", "nb positifs", "Pourcentages")

write.xlsx(final, file = "Microcin_results.xlsx", sheetName = "liste", row.names = FALSE, append = FALSE)
write.xlsx(pc, file = "Microcin_results.xlsx", sheetName = "Pourcentages", row.names = FALSE, append = TRUE)

#############################################################################################################
##################################### Analyses tous les E. coli même blast ##################################
#############################################################################################################
### Blast parameters : t = 50

### Mooving to the working directory 
setwd("~/Desktop/")

### Reading the working file
microcin = read.csv("microcin_all_fasta.csv", sep = ";", header = TRUE, stringsAsFactors = FALSE, dec = ".")
microcin$unique = gsub(".agp.fasta", "", microcin$souche)
microcin$unique = gsub(".awked.fasta", "", microcin$unique)
microcin$unique = gsub(".scfd.fasta", "", microcin$unique)

#### Number of tested strains :
length(unique(microcin$unique)) ##### 422 souches

#### Strain in which microcin have been find :
microcinpos = microcin[microcin$gene != "le nombre de genes trouves est de 0",]
microcinpos = microcinpos[!(grepl("le nombre de genes trouves est de", microcinpos$gene)),]

#### Number of unique positive strains :
length(unique(microcinpos$unique)) ##### 50 souches

#### Type of gene find :
unique(microcinpos$gene)

mcmA = microcinpos[grepl("mcmA", microcinpos$gene),]
length(unique(mcmA$unique)) #### 14 souches

mcmI = microcinpos[grepl("mcmI", microcinpos$gene),]
length(unique(mcmI$unique)) #### 14 souches

mchB = microcinpos[grepl("mchB", microcinpos$gene),]
length(unique(mchB$unique)) #### 39 souches

mchI = microcinpos[grepl("mchI", microcinpos$gene),]
length(unique(mchI$unique)) #### 39 souches

mcjA = microcinpos[grepl("mcjA", microcinpos$gene),]
length(unique(mcjA$unique)) #### 2 souches

nettoyage = function(df){
df2 = df[FALSE,]
  for (nom in unique(df$unique)){
    df1 = df[df$unique == nom,]
    sup = max(as.numeric(df1$pourcentage.d.identite))
    if (sup == 100){
      df1$pourcentage.d.identite = gsub("100.00","100", df1$pourcentage.d.identite)
    }
    print(sup)
    df1 = df1[df1$pourcentage.d.identite == sup,]
    df2 = rbind(df2,df1)
  }
  return(df2)
}

mcmA_net = nettoyage(mcmA)
mcmI_net = nettoyage(mcmI)
mchB_net = nettoyage(mchB)
mchB_net = mchB_net[!(grepl("mchB_Strain_Nissle2",mchB_net$gene)),]
mchI_net = nettoyage(mchI)
mchI_net = mchI_net[!(grepl("mchI_Strain_Nissle2",mchI_net$gene)),]



write.xlsx(unique(microcinpos$unique), file = "Microcin_results.xlsx", sheetName = "liste", row.names = FALSE, append = FALSE)
write.xlsx(mcmA_net, file = "Microcin_results.xlsx", sheetName = "mcmA", row.names = FALSE, append = TRUE)
write.xlsx(mcmI_net, file = "Microcin_results.xlsx", sheetName = "mcmI", row.names = FALSE, append = TRUE)
write.xlsx(mchB_net, file = "Microcin_results.xlsx", sheetName = "mchB", row.names = FALSE, append = TRUE)
write.xlsx(mchI_net, file = "Microcin_results.xlsx", sheetName = "mchI", row.names = FALSE, append = TRUE)
write.xlsx(mcjA, file = "Microcin_results.xlsx", sheetName = "mcjA", row.names = FALSE, append = TRUE)
