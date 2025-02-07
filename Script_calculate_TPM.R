################################################### SCRIPT CALCULO TPMs ######################################################################
setwd("C:/Users/jmper/OneDrive/Escritorio/Atlas/counts_lengths_results") # for mouse
setwd("C:/Users/jmper/OneDrive/Escritorio/Pulmones") # for human
#Librerías necesarias
library(DGEobj.utils)
library(limma)
library(edgeR)
options(scipen = 999) #Cambiar a nomenclatura no científica si Rstudio está configurado para dar los resultados en formato 10eX.
options(decimal.mark=".") #Si Rstudio no está en inglés, especificar "." como decimal y que no coja la ","


counts <- read.table(file = "sars_ricin_counts.txt",
                     sep = "\t",
                     header = TRUE,
                     row.names = 1)

counts <- as.matrix(counts)

gene_length <- read.table(file = "sars_ricin_lengts.txt",
                          sep = "\t",
                          header = TRUE,
                          row.names = 1)
gene_length <- as.matrix(gene_length)

#Calcular TPMs.
TPMs <- convertCounts(countsMatrix = counts,
                      unit = "TPM",
                      geneLength = gene_length,
                      log = FALSE,
                      normalize = "none",
                      prior.count = NULL)

#Eliminar "_X" de los nombres.
cols <- colnames(TPMs)
cols <- gsub('_.$', '', cols)
#cols <- gsub('_..$', '', cols)
cols #Comprobar que los nombres han sufrido la modificación correcta

colnames(TPMs) <- cols

#Convertir en data frame y redondear a dos dígitos.
TPMs <- as.data.frame(TPMs)
TPMs <- round(TPMs, digits = 2)

setwd("C:/Users/jmper/OneDrive/Escritorio/Atlas/counts_lengths_results")
setwd("C:/Users/jmper/OneDrive/Escritorio/Pulmones")
write.table(TPMs, "RNA-seq_SARS_ricin.txt", sep="\t", quote = FALSE, col.names=TRUE, row.names = TRUE)
  
