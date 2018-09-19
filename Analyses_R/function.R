
plasmides_trie = function(df){
  df_trie = df[FALSE,] #creation d'un tableau vide de lignes mais ayant les colones déja prêtes 
  
  for (i in 1:nrow(df)){ #on va parcourir toutes les lignes du tableau
    egalite = which(df$file == df$file[i] & 
                      df$Contig == df$Contig[i] & 
                      df$Position_in_contig == df$Position_in_contig[i])
    #Pour chaque ligne, on regarde si pour le nom de la bacterie, il y a deux contig
    #identiques qui ont ete trouves, si oui, on regarde si les positions sur le contig sont les memes
    
    if (length(egalite) <= 1){ #si deux position differentes alors length(egalite) < 1 donc 
      #on ajoute la ligne au tableau de nettoyage
      df_trie = rbind(df_trie, df[i,])
    }
    
    else{
      #si les deux positions sont identiques alors on va regarder quel est le poucentage d identite vs
      #la sequence de reference
      autre_ligne = egalite[which(egalite != i)] #contient les coordonnees des lignes identiques a le ligne i
      id_i = df$X.identity[i] #determination de la valeur du % identite pour la ligne i
      autre_id = df$X.identity[autre_ligne] #determination de la valeur du % dans les autres lignes
      
      if (id_i > max(autre_id)){
        df_trie = rbind(df_trie, df[i,]) #si le % identite est le plus fort a la ligne i alors on ajoute
        #cette ligne au tableau de nettoyage sinon on quitte
      }
    }
  }
  return(df_trie)
}

clean_antibio = function(df){
  
  df$beta.lactamase = gsub("\\([0-9]+\\.[0-9]+\\%\\,[0-9]+\\.[0-9]+\\%\\)", "", df$beta.lactamase)
  df$aminoglycoside = gsub("\\([0-9]+\\.[0-9]+\\%\\,[0-9]+\\.[0-9]+\\%\\)", "", df$aminoglycoside)
  df$phenicol = gsub("\\([0-9]+\\.[0-9]+\\%\\,[0-9]+\\.[0-9]+\\%\\)", "", df$phenicol)
  df$quinolone = gsub("\\([0-9]+\\.[0-9]+\\%\\,[0-9]+\\.[0-9]+\\%\\)", "", df$quinolone)
  df$rifampicin = gsub("\\([0-9]+\\.[0-9]+\\%\\,[0-9]+\\.[0-9]+\\%\\)", "", df$rifampicin)
  df$colistin = gsub("\\([0-9]+\\.[0-9]+\\%\\,[0-9]+\\.[0-9]+\\%\\)", "", df$colistin)
  df$fosfomycin = gsub("\\([0-9]+\\.[0-9]+\\%\\,[0-9]+\\.[0-9]+\\%\\)", "", df$fosfomycin)
  df$sulphonamide = gsub("\\([0-9]+\\.[0-9]+\\%\\,[0-9]+\\.[0-9]+\\%\\)", "", df$sulphonamide)
  df$tetracyclinet = gsub("\\([0-9]+\\.[0-9]+\\%\\,[0-9]+\\.[0-9]+\\%\\)", "", df$tetracyclinet)
  df$trimethoprim = gsub("\\([0-9]+\\.[0-9]+\\%\\,[0-9]+\\.[0-9]+\\%\\)", "", df$trimethoprim)
  df$macrolide = gsub("\\([0-9]+\\.[0-9]+\\%\\,[0-9]+\\.[0-9]+\\%\\)", "", df$macrolide)
  df$vancomycin = NULL
  
  return(df)
}

prep.itol = function(categorie, colonne){
  vec = grepl(categorie, colonne)
  vec = replace(vec, vec == "TRUE", 1)
  vec = replace(vec, vec == 0, -1)
  return(vec)
}

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


atb = function(path){
  df = read.table(path,header = TRUE, stringsAsFactors = FALSE)
  df$File = gsub("Ec", "", df$File)
  df$File = gsub("_S[0-9]+.scfd.fastq.awked", "", df$File)
  df$File = gsub(".awked", "", df$File)
  df$File = gsub("_S[0-9]+.scfd", "", df$File)
  colnames(df)
  df = df[,c("File", "beta.lactamase", "aminoglycoside", "quinolone", "fosfomycin", "trimethoprim",
             "colistin", "macrolide","phenicol","sulphonamide","fusidicacid","rifampicin","tetracyclinet",
             "nitroimidazole","oxazolidinone")]
}

mlst = function(path){
  df = read.table(path, header = TRUE, stringsAsFactors = FALSE)
  df$file = gsub("/pasteur/projets/policy01/shigella-ngs/EcCaen/Ecoli_BLSE_2018/scfd_fasta/Ec", "", df$file)
  df$file = gsub("/pasteur/projets/policy01/shigella-ngs/EcCaen/Rea_neonat/fasta/scfd_fasta/", "", df$file)
  df$file = gsub("/pasteur/projets/policy01/shigella-ngs/EcCaen/BLSE/fasta/fasta_files/awked_fasta/", "", df$file)
  df$file = gsub("_S[0-9]+.scfd.fasta", "", df$file)
  df$file = gsub(".scfd.fasta", "", df$file)
  df$file = gsub(".awked.fasta", "", df$file)
  df$allele_number = NULL
  return(df)
}

sero = function(path){
  df = read.table(path, header = TRUE, stringsAsFactors = FALSE)
  df$OH = paste(gsub("wz[a-z]_","",df$O), gsub("fliC_","",df$H), sep = ":")
  df = df[,c("file","OH")]
  df$file = gsub("/pasteur/projets/policy01/shigella-ngs/EcCaen/Ecoli_BLSE_2018/scfd_fasta/Ec", "", df$file)
  df$file = gsub("/pasteur/projets/policy01/shigella-ngs/EcCaen/Rea_neonat/fasta/awked_fasta/", "", df$file)
  df$file = gsub("/pasteur/projets/policy01/shigella-ngs/EcCaen/BLSE/fasta/fasta_files/awked_fasta/", "", df$file)
  df$file = gsub("_S[0-9]+.scfd.fasta", "", df$file)
  df$file = gsub(".awked.fasta", "", df$file)
  return(df)
}

fim = function(path){
  df = read.table(path, sep = '\t',header = TRUE, stringsAsFactors = FALSE)
  df = df[,c("Souche","Fimtype")]
  df$Souche = gsub(".agp.fasta", "", df$Souche)
  df$Souche = gsub("Ec", "", df$Souche)
  return(df)
}

viru = function(path){
df = read.table(path, header = TRUE, stringsAsFactors = FALSE)
df$virulence = apply(df[,-c(1)], 1, function(a){
  a = paste(a, collapse = ";")
  a = gsub("-;-;", "", a)
  a = gsub(";-", "", a)
  a = gsub("-;", "", a)
  a = gsub(";","-",a)
  a = gsub("-", "_",a)
  return(a)
})
df = df[,c("file","virulence")]
df$file = gsub("/pasteur/projets/policy01/shigella-ngs/EcCaen/run180223/fasta/","", df$file)
df$file = gsub("/pasteur/projets/policy01/shigella-ngs/EcCaen/Rea_neonat/fasta/awked_fasta/", "", df$file)
df$file = gsub("/pasteur/projets/policy01/shigella-ngs/EcCaen/Ecoli_BLSE_2018/scfd_fasta/Ec", "", df$file)
df$file = gsub(".scfd.fastq.awked.fasta","", df$file)
df$file = gsub(".awked.fasta","", df$file)
df$file = gsub("_S[0-9]+.scfd.fasta", "", df$file)
return(df)
}

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


## function a retravailler +++
resume_years = function(x,df){
  data = merge(resume(x),resume(x[grepl("2012",df$Date.de.prel.)]), by.x = "genes", by.y = "genes", all.x = TRUE)
  colnames(data) = c("genes", "number.total", "percentages.total", "number.2012", "percentages.2012")
  data = merge(data,resume(x[grepl("2013",df$Date.de.prel.)]), by.x = "genes", by.y = "genes", all.x = TRUE)
  colnames(data) = c("genes", "number.total", "percentages.total", "number.2012", "percentages.2012",
                     "number.2013", "percentages.2013")
  data = merge(data,resume(x[grepl("2014",df$Date.de.prel.)]), by.x = "genes", by.y = "genes", all.x = TRUE)
  colnames(data) = c("genes", "number.total", "percentages.total", "number.2012", "percentages.2012",
                     "number.2013", "percentages.2013","number.2014", "percentages.2014")
  data = merge(data,resume(x[grepl("2015",df$Date.de.prel.)]), by.x = "genes", by.y = "genes", all.x = TRUE)
  colnames(data) = c("genes", "number.total", "percentages.total","number.2012", "percentages.2012",
                     "number.2013", "percentages.2013", "number.2014", "percentages.2014","number.2015", "percentages.2015")
  data = merge(data,resume(x[grepl("2016",df$Date.de.prel.)]), by.x = "genes", by.y = "genes", all.x = TRUE)
  colnames(data) = c("genes", "number.total", "percentages.total", "number.2012", "percentages.2012",
                     "number.2013", "percentages.2013","number.2014", "percentages.2014","number.2015", 
                     "percentages.2015", "number.2016","percentages.2016")
  data = merge(data,resume(x[grepl("2017",df$Date.de.prel.)]), by.x = "genes", by.y = "genes", all.x = TRUE)
  colnames(data) = c("genes", "number.total", "percentages.total","number.2012", "percentages.2012",
                     "number.2013", "percentages.2013","number.2014", "percentages.2014","number.2015", 
                     "percentages.2015", "number.2016","percentages.2016", "number.2017", "percentages.2017")
  data[is.na(data)] = 0
  data = data[order(data$number.total, decreasing = TRUE),]
  return(data)
}