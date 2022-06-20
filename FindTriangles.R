library(data.table)
library("dplyr")
library("igraph")

#Caricamento nodi METAPATHWAY

nodes<- read.delim("metapathway_nodes.txt",
                   header = TRUE,
                   sep = "\t",
                   colClasses = "character")

#Caricamento archi METAPATHWAY

edges<- read.delim("metapathway_edges.txt",
                   header = TRUE,
                   sep = "\t",
                   colClasses = "character")

#Caricamento matrice di espressione

expressions<- read.delim("brca_expressions.txt",
                         header = TRUE,
                         sep = "\t",
                         colClasses = "character")

#Caricamento matrice di metilazione

methylations<- read.delim("brca_methylations.txt",
                          header = TRUE,
                          sep = "\t",
                          colClasses = "character")

#Caricamento matrice di mutazione

mutations<- read.delim("brca_mutations.txt",
                       header = TRUE,
                       sep = "\t",
                       colClasses = "character")


{
  
  #Sostituiamo Attivazione con 1
  edges$Subtype[edges$Subtype=="ACTIVATION"] <- 1
  
  #Sostituiamo Inibizione con -1
  edges$Subtype[edges$Subtype=="INHIBITION"] <- -1
  
  #Rinominiamo la terza colonna "weight"
  colnames(edges)[3] <- "weight"
  
  #Prendiamo in considerazione i geni che compaiono in tutti e tre i dataset
  gene_names <- intersect(intersect(intersect(nodes$Name,rownames(methylations)),expressions$Hugo_Symbol),rownames(mutations))
  
  #Rimuoviamo da tutti i dataset i record in cui compaiono i geni che non sono in comune tra tutti
  {
    expressions <- expressions[expressions$Hugo_Symbol %in% gene_names,]
    mutations <- mutations[gene_names,]
    methylations <- methylations[gene_names,]
    nodes <- nodes[nodes$Name %in% gene_names,]
    edges <- edges[edges$Start %in% nodes$Id & edges$End %in% nodes$Id,]
  }
  
  #Sistemiamo i nomi delle colonne e righe dei dataframe in modo da agevolare l'utilizzo dello stesso successivamente
  {
    rownames(expressions) <- expressions$Entrez_Gene_Id
    rownames(methylations) <- expressions$Entrez_Gene_Id
    rownames(mutations) <- expressions$Entrez_Gene_Id
    expressions$Entrez_Gene_Id <- NULL
    expressions$Hugo_Symbol <- NULL
  }
}

{
  exp_patients <- colnames(expressions)
  meth_patients <- colnames(methylations)
  mut_patients <- colnames(mutations)
  
  for(i in 1:length(exp_patients))
    exp_patients[i] <- substr(exp_patients[i],1,12)
  for(i in 1:length(meth_patients))
    meth_patients[i] <- substr(meth_patients[i],1,12)
  for(i in 1:length(mut_patients))
    mut_patients[i] <- substr(mut_patients[i],1,12)
  
  colnames(expressions) <- exp_patients
  colnames(methylations) <- meth_patients
  colnames(mutations) <- mut_patients
  
  patient_code <- intersect(intersect(exp_patients,meth_patients),mut_patients)
  
  expressions <- expressions[,patient_code]
  methylations <- methylations[,patient_code]
  mutations <- mutations[,patient_code]
}

metapathway <- graph_from_data_frame(edges,
                                     directed = TRUE,
                                     nodes)
#ggplot(metapathway)

#Creazione della matrice di adiacenza
adj_matrix <- as.matrix(as_adjacency_matrix(metapathway,
                                            names=TRUE))
n <- ncol(adj_matrix)

#Lista dei triangoli della metapathway
result <- findTriangles(adj_matrix,n)


#Data frame che conterrà per ogni paziente e per ogni triangolo, le informazioni sulla metilazione / mutazione / espressione dei geni coinvolti
output <- data.frame(matrix(NA,    # Create empty data frame
                            nrow = length(result),
                            ncol = length(colnames(expressions))))

colnames(output) <- colnames(expressions)
rownames(output) <- c(1:length(result))


#Creazione della codifica binaria di ogni triangolo e paziente
for (i in 1:length(result)){
  for(j in 1:length(colnames(expressions))){
    str <- ""
    for(k in 1:3){
      if(expressions[result[[i]][k],j]>0) str <- paste(str,"1",sep="")
      else str <- paste(str,"0",sep="")
      if(methylations[result[[i]][k],j]>0) str <- paste(str,"1",sep="")
      else str <- paste(str,"0",sep="")
      if(mutations[result[[i]][k],j]>0) str <- paste(str,"1",sep="")
      else str <- paste(str,"0",sep="")
    }
    output[i,j] <- str
  }
}

#Scrittura su disco della tabella
write.table(output,"output.txt", append = FALSE, sep = "\t",
            row.names = TRUE, col.names = TRUE)




#Funzione che restituisce una lista di vettori di tre elementi contenenti i triangoli della metapathway

findTriangles <- function(M,n) {
  out <- list()
  count <- 1
  for (i in 1:(n-1))
    for (j in (i+1):n)
      if(M[i,j]==1)
        for(k in (j+1):n)
          if(k<=n)
            if(M[j,k]==1 & M[k,i]==1 ){
              out[[count]] <- c(rownames(adj_matrix)[i],
                                rownames(adj_matrix)[j],
                                rownames(adj_matrix)[k])
              
              count <- count+1
            }
  return(out)
}

