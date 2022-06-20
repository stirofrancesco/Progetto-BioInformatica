library(data.table)
library("dplyr")
library("igraph")
library("base")
library(devtools)
library(binaryLogic)


#Caricamento matrice di codifica

data <- read.delim("output.txt",
                   header = TRUE,
                   sep = "\t",
                   colClasses = "character")

#Caricamento matrice contenenti dati sui sottotipi e stadi tumorali dei pazienti

cancer_data <- read.delim("brca_clinical.txt",
                   header = TRUE,
                   sep = "\t",
                   colClasses = "character")

{
  cancer_data <- na.omit(cancer_data)
  patients <- colnames(data)

  for (i in 1:length(patients))
    patients[i] <- gsub("[.]", "-", patients[i])

  patients <- intersect(patients,rownames(cancer_data))

  colnames(data) <- patients
  data <- data[,patients]
  cancer_data <- cancer_data[patients,]
  data<- as.data.frame(t(data))


  n_patients <- length(rownames(data))
  n_triangles <- length(colnames(data))
}


#Solo triangoli con metilazioni / mutazioni

#index <- c()
#for (i in 1:length(colnames(data))){
#  for(j in 1:length(rownames(data))){
#    print(paste(i,j))
#    if(substr(data[j,i],2,2)==0 & substr(data[j,i],5,5)==0 & substr(data[j,i],8,8)==0 & 
#       substr(data[j,i],3,3)==0 & substr(data[j,i],6,6)==0 & substr(data[j,i],9,9)==0){
#      index <- append(index,i)
#      break
#    }
#  }
#}

#data <- data[,-c(index)]


#Generazione delle configurazioni di triangoli
conf <- generate()

#Considero solo i triangoli con metilazioni/mutazioni
index <- c()
for (i in 1:length(conf)){
  if (substr(conf[i],2,2)==0 & substr(conf[i],5,5)==0 & substr(conf[i],8,8)==0 & 
      substr(conf[i],3,3)==0 & substr(conf[i],6,6)==0 & substr(conf[i],9,9)==0)
    index <- append(index,i)
}
conf <- conf[-c(index)]

#Matrice di conteggi inzialmente vuota
count <- data.frame(matrix(0,    # Create a matrix with 0s
                           nrow = n_patients,
                           ncol = length(conf)))

colnames(count) <- conf
rownames(count) <- patients

#Conteggio dei triangoli di ogni tipo
for(j in 1:length(conf)){
  count[,j] <- apply(data, 1, function(x) length(which(x==conf[j])))
  print(j) 
}

#Aggiunta delle informazioni sul sottotipo tumorale da usare per la classificazione
count$SUBTYPE <- cancer_data$SUBTYPE

#Scrittura della matrice dei conteggi su disco
write.table(count,"count.txt", append = FALSE, sep = "\t",
            row.names = TRUE, col.names = TRUE)


#Funzione che genera tutte le combinazioni possibili di triangoli
generate <- function(){
  config <- c()
  for (i in 1:512) {
    config[i] <- substr(paste(rev(as.integer(intToBits(i-1))), collapse=""),24,32)
  }
  config 
}
