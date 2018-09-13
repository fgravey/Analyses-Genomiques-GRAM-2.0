## Function definitions
plasmides_asso = function(x){
  return(paste(sort(plasmide$Plasmide[grepl(x, plasmide$Strain)]), collapse = "_"))
}

### Reading files 
epi = read.csv("/Volumes/Maxtor/Back_up_pasteur/Kp_Caen/data/kp_epi.csv", header = TRUE, sep = ";",
               stringsAsFactors = FALSE)
atb = read.csv("/Volumes/Maxtor/Back_up_pasteur/Kp_Caen/Analyses/Resistance/liste.Resfinder.csv", header = TRUE,
               sep = ';', stringsAsFactors = FALSE)
mlst = read.csv("/Volumes/Maxtor/Back_up_pasteur/Kp_Caen/Analyses/mlst/mlst.csv", 
                header = TRUE, sep = ";",stringsAsFactors = FALSE)

plasmide = read.csv("/Volumes/Maxtor/Back_up_pasteur/Kp_Caen/Analyses/plasmides_2/plasmidfinder_results.csv",
                    header = TRUE, sep = ";", stringsAsFactors = FALSE)

### Keeping only relevant informations
epi = epi[,-c(10:66)]
atb = atb[,-c(4,6,8,9,12,16)]

##merging the df
kp = merge(epi, mlst, by.x = "N_etude", by.y = "N_etude", all.y = TRUE)
kp = merge(kp,atb, by.x = "N_etude", by.y = "File")

###creation of the plasmides column
plasm = as.vector(sapply(sort(unique(plasmide$Strain)), plasmides_asso))
kp$plasmides = plasm

## wrtiting the output file
write.table(kp,"/Volumes/Maxtor/Back_up_pasteur/Kp_Caen/Analyses/kp_avant_2018.csv", sep = ";", row.names = FALSE)


