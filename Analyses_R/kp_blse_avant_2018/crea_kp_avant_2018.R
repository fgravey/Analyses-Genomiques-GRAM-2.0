### Reading files 
epi = read.csv("/Volumes/Maxtor/Back_up_pasteur/Kp_Caen/data/kp_epi.csv", header = TRUE, sep = ";",
               stringsAsFactors = FALSE)
atb = read.csv("/Volumes/Maxtor/Back_up_pasteur/Kp_Caen/Analyses/Resistance/liste.Resfinder.csv", header = TRUE,
               sep = ';', stringsAsFactors = FALSE)

### Keeping only relevant informations
epi = epi[,-c(10:66)]
atb = atb[,-c(4,6,8,9,12,16)]

##merging the df
kp = merge(epi, atb, by.x = "N_etude", by.y = "File", all.y = TRUE)

## wrtiting the output file
write.table(kp,"/Volumes/Maxtor/Back_up_pasteur/Kp_Caen/Analyses/kp_avant_2018.csv", sep = ";", row.names = FALSE)
