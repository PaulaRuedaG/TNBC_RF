library(Biobase)
library(limma)    

datos <- read.csv("path/fRMA_normalized_expression.csv", row.names = 1)  # Usa row.names = 1 si los genes están en la primera columna

colnames(datos) <- gsub("_.*", "", colnames(datos)) 

muestras_TNBC_RD1 <- readLines("path/muestras_TNBC_RD_GSE25066.txt") 
muestras_TNBC_pCR1 <- readLines("path/muestras_TNBC_pCR_GSE25066.txt") 
muestras_TNBC_RD2 <- readLines("path/muestras_TNBC_RD_GSE20194.txt") 
muestras_TNBC_pCR2 <- readLines("path/muestras_TNBC_pCR_GSE20194.txt") 
muestras_TNBC_RD <- c(muestras_TNBC_RD1,muestras_TNBC_RD2)
muestras_TNBC_pCR <- c(muestras_TNBC_pCR1,muestras_TNBC_pCR2)
muestras_seleccionadas <- c(muestras_TNBC_RD,muestras_TNBC_pCR)


grupo <- ifelse(colnames(datos) %in% muestras_TNBC_RD, "RD", 
                ifelse(colnames(datos) %in% muestras_TNBC_pCR, "pCR", NA))  # NA si alguna no coincide

metadata <- data.frame(Sample = colnames(datos), Group = grupo)
metadata <- metadata[!is.na(metadata$Group), ] 
metadata <- metadata[metadata$Sample %in% colnames(datos), ]
datos <- datos[, colnames(datos) %in% metadata$Sample]

metadata <- metadata[match(colnames(datos), metadata$Sample), ]

stopifnot(all(colnames(datos) == metadata$Sample))

rownames(metadata) <- metadata$Sample  
metadata$Sample <- NULL  
metadata <- new("AnnotatedDataFrame", data = metadata)  

assay_data <- assayDataNew("environment", exprs = as.matrix(datos))

expSet <- ExpressionSet(assayData = as.matrix(datos), phenoData = metadata)

expSet
pData(expSet)
exprs(expSet)

realizar_bootstrap_DEGs <- function(muestras_TNBC_RD, muestras_TNBC_pCR, expSet, num_iter = num_iteraciones, logFC_threshold, p_value_threshold) {
  
  lista_acumulada_DEGs <- list()
  
  for (i in 1:num_iter) {
    cat("Iteración: ", i, "\n")
    
    muestras_RD_bootstrap <- sample(muestras_TNBC_RD2, 50, replace = TRUE)
    muestras_pCR_bootstrap <- sample(muestras_TNBC_pCR2, 50, replace = TRUE)
    muestras_bootstrap <- c(muestras_RD_bootstrap, muestras_pCR_bootstrap) 
    datos_bootstrap <- expSet[, muestras_bootstrap]
    
    grupos_bootstrap <- factor(c(rep("RD", length(muestras_RD_bootstrap)), rep("pCR", length(muestras_pCR_bootstrap))))
    
    design_bootstrap <- model.matrix(~0 + grupos_bootstrap)
    colnames(design_bootstrap) <- c("RD", "pCR")
    
    fit_bootstrap <- lmFit(exprs(datos_bootstrap), design_bootstrap)
    contraste_bootstrap <- makeContrasts(RD_vs_pCR = RD - pCR, levels = design_bootstrap)
    fit_bootstrap2 <- contrasts.fit(fit_bootstrap, contrasts = contraste_bootstrap)
  
    fit_bootstrap2 <- eBayes(fit_bootstrap2)
    
    tT_bootstrap <- topTable(fit_bootstrap2, coef = 1, adjust = "fdr", number = Inf)
    
    DEGs_filtrados <- tT_bootstrap[abs(tT_bootstrap$logFC) > logFC_threshold & tT_bootstrap$adj.P.Val < p_value_threshold, ]
    
    lista_acumulada_DEGs[[i]] <- rownames(DEGs_filtrados)
  }
  
  DEGs_combinados <- unlist(lista_acumulada_DEGs)
  
  return(DEGs_combinados)
}

#---------------------------------------------------
# Bootstrap
#---------------------------------------------------

num_iteraciones <- 999  # Número de iteraciones del bootstrap
logFC_threshold <- 1  # Umbral para logFC
p_value_threshold <- 0.05  # Umbral de significancia estadística

# Ejecutar la función y almacenar todos los resultados
resultados_bootstrap <- realizar_bootstrap_DEGs(muestras_TNBC_RD, muestras_TNBC_pCR, expSet, num_iteraciones, logFC_threshold, p_value_threshold)

#---------------------------------------------------
# DEGs storage
#---------------------------------------------------

DEGs_combinados <- unlist(resultados_bootstrap)
tabla_frecuencias_DEGs <- as.data.frame(table(DEGs_combinados))
tabla_frecuencias_DEGs <- tabla_frecuencias_DEGs[order(-tabla_frecuencias_DEGs$Freq), ]
head(tabla_frecuencias_DEGs)

# Seleccionar los DEGs presentes en n% de las iteraciones
umbral_frecuencia <- 0.2  # Umbral de frecuencia (30%)
genes_destacados <- tabla_frecuencias_DEGs$DEGs_combinados[tabla_frecuencias_DEGs$Freq >= (num_iteraciones * umbral_frecuencia)]

#---------------------------------------------------
# Volcano plot 
#---------------------------------------------------

grupos <- factor(c(rep("RD", length(muestras_TNBC_RD)), rep("pCR", length(muestras_TNBC_pCR))))
design <- model.matrix(~0 + grupos)
colnames(design) <- c("RD", "pCR")

fit <- lmFit(exprs(expSet), design)
contraste <- makeContrasts(RD_vs_pCR = RD - pCR, levels = design)
fit2 <- contrasts.fit(fit, contrasts = contraste)
fit2 <- eBayes(fit2)

resultados <- topTable(fit2, coef=1, number=Inf, adjust="fdr")

ids_sondas <- rownames(resultados)

genes_destacados_character <- as.character(genes_destacados)
highlight <- which(ids_sondas %in% genes_destacados_character)

logFC_values <- resultados$logFC
p_values <- resultados$P.Value
fdr_values <- resultados$adj.P.Val

volcanoplot(fit2, coef=1, main="Volcano Plot RD vs pCR", pch=20,
            highlight=highlight, col=ifelse(1:nrow(fit2) %in% highlight, "red", "black"),
            xlab = "Log2 Fold Change", ylab = "-Log10(Adjusted P-value)")

highlight_logFC <- logFC_values[highlight]
highlight_p_value <- -log10(p_values[highlight])

text(highlight_logFC, highlight_p_value, labels = ids_sondas[highlight], 
     pos = 3, col = "red", cex = 0.7)
legend("topright", legend = c("Probes destacados"), col = "red", pch = 20)

cat("IDs de las sondas resaltadas en el volcano plot:\n")
candidatos <- ids_sondas[highlight]
print(candidatos)

relacion_sonda_gen <- read.csv("path/probes_genes.csv")

tabla_exportar <- data.frame(
  ID = ids_sondas,
  GeneSymbol = relacion_sonda_gen$GeneSymbol[match(ids_sondas, relacion_sonda_gen$ID)],
  LogFC = logFC_values,
  p_values = p_values,
  fdr_values = fdr_values
)

print(tabla_exportar)

write.csv(tabla_exportar, file = "genes_destacados_volcano.csv", row.names = FALSE)

mtrx_expresion=exprs(expSet)
mtrx_expresion_filtrada <- as.data.frame(t(mtrx_expresion[, muestras_seleccionadas, drop = FALSE]))
mtrx_expresion_filtrada$grupos <- grupos
degs_filtrados<- mtrx_expresion_filtrada[, c(genes_interes, "grupos"), drop = FALSE]
write.csv(degs_filtrados, "datos_ML_Conjuntos.csv", row.names = TRUE)
