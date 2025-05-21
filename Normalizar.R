if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "affy", "oligo","hgu133plus2frmavecs","frma"))

library(GEOquery)
library(affy)   # Para procesar Affymetrix CEL files
library(oligo)  # Para arrays Affymetrix de nueva generación
library(frma)
library(hgu133plus2frmavecs)  # Si usas Affymetrix HGU133 Plus 2.0
library(hgu133afrmavecs)

gse_ids <- c("GSE25066", "GSE20194")

for (gse in gse_ids) {
  getGEOSuppFiles(gse)  # Descarga archivos suplementarios (incluye los CEL)
  untar(paste0(gse, "/", gse, "_RAW.tar"), exdir = paste0(gse, "_CEL"))  # Extrae los CEL
}

data1 <- ReadAffy(celfile.path = "C:/Users/pasan/Documents/Maestria/Trabajo_de_grado/GSE25066/")
data2 <- ReadAffy(celfile.path = "C:/Users/pasan/Documents/Maestria/Trabajo_de_grado/GSE20194/")

frma_data1 <- frma(data1)
frma_data2 <- frma(data2)

exprs1 <- exprs(frma_data1)
exprs2 <- exprs(frma_data2)

common_genes <- intersect(rownames(exprs1), rownames(exprs2))

exprs1 <- exprs1[common_genes, ]
exprs2 <- exprs2[common_genes, ]

combined_exprs <- cbind(exprs1, exprs2)

pheno_data <- rbind(pData(frma_data1), pData(frma_data2)) 
combined_eset <- ExpressionSet(assayData = combined_exprs, phenoData = AnnotatedDataFrame(pheno_data))

write.csv(exprs(combined_eset), "fRMA_normalized_expression.csv")

