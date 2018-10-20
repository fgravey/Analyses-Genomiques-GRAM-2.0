## Working library
library(lubridate)
library(ggplot2)

### function definition

compteur = function(df,annee_1,annee_2){
  date = c()
  nb = c()
for (i in seq(annee_1,annee_2)){
  for (j in seq(1,12)){
    annee= df[year(df$Date.prelevement) == i,]
    mois = sum(month(annee$Date.prelevement) == j)
    print(mois)
    nb = append(nb, mois)
    if (j < 10){
      j = paste(0,j,sep = "")
    }
    date = append(date, paste(i, j,"01", sep = "/"))
  }
}

out = data.frame(Dates = ymd(date), Nombres.kpblse= nb, stringsAsFactors = FALSE)

return(out)
}

#### Mooving to the working directory
setwd("/Volumes/Maxtor/Back_up_pasteur/kp_caen_2018/Data/epi/")

## Reading the csv file and cleaning data
blse1218 = read.csv("blse_caen_1218.csv",header = TRUE, stringsAsFactors = FALSE, sep = ";")

##############################################################################################################
################################ Working on kpn #################################################
##############################################################################################################
kp1218 = blse1218[blse1218$Germe == "Klebsiella pneumoniae pneumoniae",]

kp1218 = kp1218[kp1218$Correspondant != "CH JACQUES MONOD" &
                  kp1218$Correspondant != "SERVICE :" & kp1218$Correspondant != "LAM BIONACRE -  VAUCELLES"
                & kp1218$Correspondant != "CSR KORIAN COTE NORMANDE" & 
                  kp1218$Correspondant !="LAM DES ANDAINES-NORMABIO",]

kp1218$Date.prelevement = gsub("/", "", kp1218$Date.prelevement)
kp1218$Date.prelevement = dmy(kp1218$Date.prelevement)

## Keeping non doublon patient per year

kp12 = kp1218[year(kp1218$Date.prelevement) == 2012,]
kp12 = kp12[!duplicated(kp12$Patient),]

kp13 = kp1218[year(kp1218$Date.prelevement) == 2013,]
kp13 = kp13[!duplicated(kp13$Patient),]

kp14 = kp1218[year(kp1218$Date.prelevement) == 2014,]
kp14 = kp14[!duplicated(kp14$Patient),]

kp15 = kp1218[year(kp1218$Date.prelevement) == 2015,]
kp15 = kp15[!duplicated(kp15$Patient),]

kp16 = kp1218[year(kp1218$Date.prelevement) == 2016,]
kp16 = kp16[!duplicated(kp16$Patient),]

kp17 = kp1218[year(kp1218$Date.prelevement) == 2017,]
kp17 = kp17[!duplicated(kp17$Patient),]

kp18 = kp1218[year(kp1218$Date.prelevement) == 2018,]
kp18 = kp18[!duplicated(kp18$Patient),]

### merging all the df
### Important !!!! Automatiser la fonction rbind !!! 
x.n = c("kp12","kp13", "kp14", "kp15", "kp16", "kp17", "kp18")
x.list = lapply(x.n, get)
kpblse1218 = do.call(rbind,x.list)

################################# Compte du nombre de kpblse par mois grâce à la fonction compteur ########
nb.kpblse.1218 = compteur(kpblse1218,2012,2018)
colnames(nb.kpblse.1218) = c("Dates", "Nombres.kpblse")
sum(nb.kpblse.1218$Nombres.kpblse)

##############################################################################################################
################################ Working on kpn #################################################
##############################################################################################################
ecol = blse1218[blse1218$Germe == "Escherichia coli",]
ecol = ecol[grepl("^[0-9]+$", ecol$Correspondant),]

ecol$Date.prelevement = gsub("/", "", ecol$Date.prelevement)
ecol$Date.prelevement = dmy(ecol$Date.prelevement)

## Keeping non doublon patient per year

ecol12 = ecol[year(ecol$Date.prelevement) == 2012,]
ecol12 = ecol12[!duplicated(ecol12$Patient),]

ecol13 = ecol[year(ecol$Date.prelevement) == 2013,]
ecol13 = ecol13[!duplicated(ecol13$Patient),]

ecol14 = ecol[year(ecol$Date.prelevement) == 2014,]
ecol14 = ecol14[!duplicated(ecol14$Patient),]

ecol15 = ecol[year(ecol$Date.prelevement) == 2015,]
ecol15 = ecol15[!duplicated(ecol15$Patient),]

ecol16 = ecol[year(ecol$Date.prelevement) == 2016,]
ecol16 = ecol16[!duplicated(ecol16$Patient),]

ecol17 = ecol[year(ecol$Date.prelevement) == 2017,]
ecol17 = ecol17[!duplicated(ecol17$Patient),]

ecol18 = ecol[year(ecol$Date.prelevement) == 2018,]
ecol18 = ecol18[!duplicated(ecol18$Patient),]


### rbind all the data frames
x.n = c("ecol12","ecol13", "ecol14", "ecol15", "ecol16", "ecol17", "ecol18")
x.list = lapply(x.n, get)
ecolblse1218 = do.call(rbind,x.list)

################################# Compte du nombre de ecolblse par mois grâce à la fonction compteur ########
nb.ecolblse.1218 = compteur(ecolblse1218,2012,2018)
colnames(nb.ecolblse.1218) = c("Dates", "Nombres.ecolblse")
sum(nb.ecolblse.1218$Nombres.ecolblse)


##############################################################################################################
################################ Working on demande per year #################################################
##############################################################################################################
## Reading the csv file and cleaning data
demandes1218 = read.csv("demandes_blse_caen_1218.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
demandes1218 = demandes1218[grepl("^[0-9]+$", demandes1218$Correspondant),]
demandes1218$Date.prelevement = dmy(demandes1218$Date.prelevement)
demandes1218 = demandes1218[order(demandes1218$Date.prelevement),]

## Keeping non doublon patient per year
demandes12 = demandes1218[year(demandes1218$Date.prelevement) == 2012,]
demandes12 = demandes12[!duplicated(demandes12$Patient),]

demandes13 = demandes1218[year(demandes1218$Date.prelevement) == 2013,]
demandes13 = demandes13[!duplicated(demandes13$Patient),]

demandes14 = demandes1218[year(demandes1218$Date.prelevement) == 2014,]
demandes14 = demandes14[!duplicated(demandes14$Patient),]

demandes15 = demandes1218[year(demandes1218$Date.prelevement) == 2015,]
demandes15 = demandes15[!duplicated(demandes15$Patient),]

demandes16 = demandes1218[year(demandes1218$Date.prelevement) == 2016,]
demandes16 = demandes16[!duplicated(demandes16$Patient),]

demandes17 = demandes1218[year(demandes1218$Date.prelevement) == 2017,]
demandes17 = demandes17[!duplicated(demandes17$Patient),]

demandes18 = demandes1218[year(demandes1218$Date.prelevement) == 2018,]
demandes18 = demandes18[!duplicated(demandes18$Patient),]

### merging all the df
demandes = c("demandes12","demandes13","demandes14","demandes15","demandes16","demandes17","demandes18")
demandes.list = lapply(demandes, get)

########################################### demandes.sans.doublon contient les patients par an sans doublon ##
demandes.sans.doublon = do.call(rbind, demandes.list)

################################# Compte du nombre de demandes par mois grâce à la fonction compteur ########
nb.demandes.1218 = compteur(demandes.sans.doublon,2012,2018)
colnames(nb.demandes.1218) = c("Dates", "Nombres.demandes")

##############################################################################################################
################################ Summary of the data ########################################################
##############################################################################################################

nb.demandes.rea = demandes.sans.doublon[grepl("^1610$|^1630$|^611$|^612$|^620$",demandes.sans.doublon$Correspondant),]
ecolrea = ecolblse1218[grepl("^1610$|^1630$|^611$|^612$|^620$",ecolblse1218$Correspondant),]
kprea = kpblse1218[grepl("^1610$|^1630$|^611$|^612$|^620$",kpblse1218$Correspondant),]

nbecolrea = compteur(ecolrea,2012,2018)
colnames(nbecolrea) = c("Dates", "Nombres.ecolblse")
nbkprea = compteur(kprea,2012,2018)
colnames(nbkprea) = c("Dates", "Nombres.kpblse")
rea1218 = compteur(nb.demandes.rea,2012,2018)
colnames(rea1218) = c("Dates", "Nombres.demandes")

rea = merge(rea1218, nbecolrea, by.x = "Dates", by.y = "Dates")
rea$tauxecol = round(rea$Nombres.ecolblse/rea$Nombres.demandes, digits = 2)
rea = merge(rea, nbkprea, by.x = "Dates", by.y = "Dates")
rea$tauxkp = round(rea$Nombres.kpblse/rea$Nombres.demandes, digits = 2)
rea[is.na(rea)] = 0


setwd("~/Desktop/")

resume = merge(nb.demandes.1218, nb.kpblse.1218, by.x = "Dates", by.y = "Dates")
resume$tauxkpn = round(resume$Nombres.kpblse/resume$Nombres.demandes, digits = 2)
resume = merge(resume, nb.ecolblse.1218, by.x = "Dates", by.y = "Dates")
resume$tauxecol = round(resume$Nombres.ecolblse/resume$Nombres.demandes, digits = 2)
resume[is.na(resume)] = 0
write.table(ecolrea, "ecolrea.csv", sep = ";", row.names = FALSE)
write.table(kprea, "kprea.csv", sep = ";", row.names = FALSE)
write.table(nb.demandes.rea, "demandes_rea.csv", sep = ";", row.names = FALSE)
write.table(rea, "resume_rea.csv", sep = ";", row.names = FALSE)



ggplot(data=rea) +
  geom_line(aes(x=Dates, y = Nombres.kpblse, color = "blse"))+
  geom_line(aes(x=Dates, y = Nombres.ecolblse, color = "ecol"))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y", limits = c(min(resume$Dates), max(resume$Dates)))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top")

ggplot(data=rea) +
  #geom_line(aes(x=Dates, y = Nombres.kpblse, color = "kpblse"))+
  geom_line(aes(x=Dates, y = tauxkp, color = "taux kpn"))+
  geom_line(aes(x=Dates, y = tauxecol, color = "taux ecol"))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y", limits = c(min(resume$Dates), max(resume$Dates)))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top")


plot(resume$Dates,resume$Nombres.demandes, type = "l")

ggplot(data=resume) +
  #geom_line(aes(x=Dates, y = Nombres.kpblse, color = "kpblse"))+
  geom_line(aes(x=Dates, y = tauxkpn, color = "taux kpn"))+
  geom_line(aes(x=Dates, y = tauxecol, color = "taux ecol"))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y", limits = c(min(resume$Dates), max(resume$Dates)))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top")

ggplot(data=resume) +
  geom_line(aes(x=Dates, y = Nombres.kpblse, color = "blse"))+
  geom_line(aes(x=Dates, y = Nombres.ecolblse, color = "ecol"))+
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y", limits = c(min(resume$Dates), max(resume$Dates)))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "top")

sum(year(demandes1218$Date.prelevement == 2014))
