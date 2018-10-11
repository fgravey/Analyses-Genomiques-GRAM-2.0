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
              