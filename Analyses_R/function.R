
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
  df$file = gsub("_S[0-9]+.scfd.fasta", "", df$file)
  df$file = gsub(".scfd.fasta", "", df$file)
  df$allele_number = NULL
  return(df)
}