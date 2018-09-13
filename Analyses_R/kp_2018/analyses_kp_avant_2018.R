## libraries
library(xlsx)

## function definition
clean_percent = function(x){
  return(gsub("\\([0-9]+\\.[0-9]+\\%\\,[0-9]+\\.[0-9]+\\%\\)","",x))
}

distribution = function(x){
  var = paste(x, collapse = "_")
  var = strsplit(var, "_")[[1]]
  return(sort(table(var), decreasing = TRUE))
}

resume = function(x){
  df = data.frame(names(distribution(x)),as.vector(distribution(x)), round(as.vector(prop.table(distribution(x)))*100,digit =2))
  colnames(df) = c("genes", "number", "percentages")
  return(df)
}

## Reading the information file
kp = read.csv("/Volumes/Maxtor/Back_up_pasteur/Kp_Caen/Analyses/kp_avant_2018.csv",
              sep = ";", header = TRUE, stringsAsFactors = FALSE)

## Cleaning all the percentage
kp = as.data.frame(t(apply(kp, 1, clean_percent)))

##### Writing files
## moving to the working directory
setwd("/Volumes/Maxtor/Back_up_pasteur/Kp_Caen/Analyses/stats")

write.xlsx(resume(kp$MLST), file = "stats_kp_avant_2018.xlsx",
           sheetName = "mlst", row.names = FALSE, append = FALSE)
write.xlsx(resume(kp$plasmides), file = "stats_kp_avant_2018.xlsx",
           sheetName = "plasmides", row.names = FALSE, append = TRUE)
write.xlsx(resume(kp$aminoglycoside), file = "stats_kp_avant_2018.xlsx",
           sheetName = "aminosides", row.names = FALSE, append = TRUE)
write.xlsx(resume(kp$beta.lactamase), file = "stats_kp_avant_2018.xlsx",
           sheetName = "beta_lactamases", row.names = FALSE, append = TRUE)
write.xlsx(resume(kp$quinolone), file = "stats_kp_avant_2018.xlsx",
           sheetName = "quinolones", row.names = FALSE, append = TRUE)
write.xlsx(resume(kp$fosfomycin), file = "stats_kp_avant_2018.xlsx",
           sheetName = "fosfomycine", row.names = FALSE, append = TRUE)
write.xlsx(resume(kp$trimethoprim), file = "stats_kp_avant_2018.xlsx",
           sheetName = "trimethoprim", row.names = FALSE, append = TRUE)
write.xlsx(resume(kp$macrolide), file = "stats_kp_avant_2018.xlsx",
           sheetName = "macrolides", row.names = FALSE, append = TRUE)
write.xlsx(resume(kp$phenicol), file = "stats_kp_avant_2018.xlsx",
           sheetName = "phenicol", row.names = FALSE, append = TRUE)
write.xlsx(resume(kp$sulphonamide), file = "stats_kp_avant_2018.xlsx",
           sheetName = "sulphonamide", row.names = FALSE, append = TRUE)
write.xlsx(resume(kp$tetracyclinet), file = "stats_kp_avant_2018.xlsx",
           sheetName = "tetracycline", row.names = FALSE, append = TRUE)

resume(kp$MLST)
resume(kp$MLST[grepl("2013",kp$Date.de.prel.)])
resume(kp$MLST[grepl("2014",kp$Date.de.prel.)])
resume(kp$MLST[grepl("2015",kp$Date.de.prel.)])
resume(kp$MLST[grepl("2016",kp$Date.de.prel.)])
resume(kp$MLST[grepl("2017",kp$Date.de.prel.)])

merge(resume(kp$MLST),resume(kp$MLST[grepl("2013",kp$Date.de.prel.)]), by.x = "genes", by.y = "genes", all.x = TRUE)
resume_years = function(x){
  data = merge(resume(x),resume(x[grepl("2013",kp$Date.de.prel.)]), by.x = "genes", by.y = "genes", all.x = TRUE)
  colnames(data) = c("genes", "number.total", "percentages.total", "number.2013", "percentages.2013")
  data = merge(data,resume(x[grepl("2014",kp$Date.de.prel.)]), by.x = "genes", by.y = "genes", all.x = TRUE)
  colnames(data) = c("genes", "number.total", "percentages.total", "number.2013", "percentages.2013",
                     "number.2014", "percentages.2014")
  data = merge(data,resume(x[grepl("2015",kp$Date.de.prel.)]), by.x = "genes", by.y = "genes", all.x = TRUE)
  colnames(data) = c("genes", "number.total", "percentages.total", "number.2013", "percentages.2013",
                     "number.2014", "percentages.2014","number.2015", "percentages.2015")
  data = merge(data,resume(x[grepl("2016",kp$Date.de.prel.)]), by.x = "genes", by.y = "genes", all.x = TRUE)
  colnames(data) = c("genes", "number.total", "percentages.total", "number.2013", "percentages.2013",
                     "number.2014", "percentages.2014","number.2015", "percentages.2015", "number.2016",
                     "percentages.2016")
  data = merge(data,resume(x[grepl("2017",kp$Date.de.prel.)]), by.x = "genes", by.y = "genes", all.x = TRUE)
  colnames(data) = c("genes", "number.total", "percentages.total", "number.2013", "percentages.2013",
                     "number.2014", "percentages.2014","number.2015", "percentages.2015", "number.2016",
                     "percentages.2016", "number.2017", "percentages.2017")
  data[is.na(data)] = 0
  data = data[order(data$number.total, decreasing = TRUE),]
  return(data)
}

write.xlsx(resume_years(kp$MLST), file = "stats_kp_avant_2018_essai_2.xlsx",
           sheetName = "mlst", row.names = FALSE, append = FALSE)
write.xlsx(resume_years(kp$plasmides), file = "stats_kp_avant_2018_essai_2.xlsx",
           sheetName = "plasmides", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(kp$aminoglycoside), file = "stats_kp_avant_2018_essai_2.xlsx",
           sheetName = "aminosides", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(kp$beta.lactamase), file = "stats_kp_avant_2018_essai_2.xlsx",
           sheetName = "beta_lactamases", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(kp$quinolone), file = "stats_kp_avant_2018_essai_2.xlsx",
           sheetName = "quinolones", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(kp$fosfomycin), file = "stats_kp_avant_2018_essai_2.xlsx",
           sheetName = "fosfomycine", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(kp$trimethoprim), file = "stats_kp_avant_2018_essai_2.xlsx",
           sheetName = "trimethoprim", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(kp$macrolide), file = "stats_kp_avant_2018_essai_2.xlsx",
           sheetName = "macrolides", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(kp$phenicol), file = "stats_kp_avant_2018_essai_2.xlsx",
           sheetName = "phenicol", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(kp$sulphonamide), file = "stats_kp_avant_2018_essai_2.xlsx",
           sheetName = "sulphonamide", row.names = FALSE, append = TRUE)
write.xlsx(resume_years(kp$tetracyclinet), file = "stats_kp_avant_2018_essai_2.xlsx",
           sheetName = "tetracycline", row.names = FALSE, append = TRUE)

essai = resume_years(kp$fosfomycin)


