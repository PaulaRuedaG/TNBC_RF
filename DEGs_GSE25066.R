#----------------------------------------------------------
#Package install 
#----------------------------------------------------------

install_and_load <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)  # Instalar paquete si no está instalado
    library(package, character.only = TRUE)  # Cargar paquete
  }
}

bioc_install_and_load <- function(package) {
  if (!require(package, character.only = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")  # Instalar BiocManager si no está instalado
    BiocManager::install(package)  # Instalar paquete desde Bioconductor
    library(package, character.only = TRUE)  # Cargar paquete
  }
}

bioc_install_and_load("GEOquery")       # Paquete para obtener datos de GEO
bioc_install_and_load("limma")          # Paquete para análisis de datos de microarrays
bioc_install_and_load("rafalib")        # Paquete para funciones adicionales de gráficos y análisis
bioc_install_and_load("EnhancedVolcano")# Paquete para plotear DEGs (genes diferencialmente expresados)
install_and_load("gplots")              # Paquete para funciones de gráficos adicionales
install_and_load("pheatmap")            # Paquete para generar heatmaps
install_and_load("ggplot2")             # Paquete para visualización de datos

#----------------------------------------------------------
#Download and pre-processing 
#----------------------------------------------------------

datos_GSE25066 <- getGEO("GSE25066", GSEMatrix = TRUE)
ID_plataforma <- levels(as.factor((datos_GSE25066[[1]])$platform_id))
if (length(datos_GSE25066) > 1) {
  idx <- grep("GPL96", attr(datos_GSE25066, "names"))
} else {
  idx <- 1
}
datos_GSE25066 <- datos_GSE25066[[idx]]
expresion_genica <- exprs(datos_GSE25066)
quantiles <- as.numeric(quantile(expresion_genica, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))

# Logarithmic transformation 
ex <- exprs(datos_GSE25066)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) #Si la diferencia entre el valor máximo (100%) y el mínimo (0%) es mayor a 50 y el primer cuartil es mayor que 0, los datos tienen una dispersión muy amplia
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(datos_GSE25066) <- log2(ex) }

#----------------------------------------------------------
#Group prepare for DEGs analysis 
#----------------------------------------------------------

# Get features to filter 
caracteristicas <- pData(datos_GSE25066)
info <- as.data.frame(caracteristicas)
info <- info[, -c(70:81)]
rows_to_modify <- 311:508
#
info[rows_to_modify, 23] <- info[rows_to_modify, 22]
info[rows_to_modify, 22] <- info[rows_to_modify, 21]
info[rows_to_modify, 21] <- info[rows_to_modify, 20]
info[rows_to_modify, 20] <- info[rows_to_modify, 19]
info[rows_to_modify, 19] <- info[rows_to_modify, 18]
info[rows_to_modify, 18] <- info[rows_to_modify, 17]
info[rows_to_modify, 17] <- info[rows_to_modify, 16]
info[rows_to_modify, (16)] <- ""

columns_to_clean <- c("characteristics_ch1.4", "characteristics_ch1.5", "characteristics_ch1.3","characteristics_ch1.11","characteristics_ch1.18")
info[columns_to_clean] <- lapply(info[columns_to_clean], function(x) gsub(".*: ", "", x))

#Select samples for RD and pCR 
TNBC_RD <- subset(info, characteristics_ch1.3 == "N" & characteristics_ch1.4 == "N" & characteristics_ch1.5 == "N" & characteristics_ch1.11 == "RD")
TNBC_pCR <- subset(info, characteristics_ch1.3 == "N" & characteristics_ch1.4 == "N" & characteristics_ch1.5 == "N" & characteristics_ch1.11 == "pCR")

muestras_TNBC_RD <- rownames(TNBC_RD)
muestras_TNBC_pCR <- rownames(TNBC_pCR)
muestras_seleccionadas <- c(muestras_TNBC_RD, muestras_TNBC_pCR)
datos_TNBC <- datos_GSE25066[, muestras_seleccionadas]

#-------------------------------------------------------
# Bootstrap function 
#-------------------------------------------------------
realizar_bootstrap_DEGs <- function(muestras_TNBC_RD, muestras_TNBC_pCR, datos_TNBC, num_iter = num_iteraciones, logFC_threshold = 1, p_value_threshold = 0.05) {
  
  lista_acumulada_DEGs <- list()
  
  for (i in 1:num_iter) {
    cat("Iteración: ", i, "\n")
    
    muestras_RD_bootstrap <- sample(muestras_TNBC_RD, 50, replace = TRUE)
    muestras_pCR_bootstrap <- sample(muestras_TNBC_pCR, 50, replace = TRUE)
    muestras_bootstrap <- c(muestras_RD_bootstrap, muestras_pCR_bootstrap) 
    datos_bootstrap <- datos_TNBC[, muestras_bootstrap]
    
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
# Bootstrap execution 
#---------------------------------------------------

num_iteraciones <- 999  # Número de iteraciones del bootstrap
logFC_threshold <- 1  # Umbral para logFC
p_value_threshold <- 0.05  # Umbral de significancia estadística

# Ejecutar la función y almacenar todos los resultados
resultados_bootstrap <- realizar_bootstrap_DEGs(muestras_TNBC_RD, muestras_TNBC_pCR, datos_TNBC, num_iteraciones, logFC_threshold, p_value_threshold)

#---------------------------------------------------
# DEGs storage
#---------------------------------------------------

DEGs_combinados <- unlist(resultados_bootstrap)
tabla_frecuencias_DEGs <- as.data.frame(table(DEGs_combinados))
tabla_frecuencias_DEGs <- tabla_frecuencias_DEGs[order(-tabla_frecuencias_DEGs$Freq), ]
head(tabla_frecuencias_DEGs)

# Seleccionar los DEGs presentes en n% de las iteraciones
umbral_frecuencia <- 0.2  # Umbral de frecuencia (10%)
genes_destacados <- tabla_frecuencias_DEGs$DEGs_combinados[tabla_frecuencias_DEGs$Freq >= (num_iteraciones * umbral_frecuencia)]

#---------------------------------------------------
# Volcano plot
#---------------------------------------------------

# Experimental design 
grupos <- factor(c(rep("RD", length(muestras_TNBC_RD)), rep("pCR", length(muestras_TNBC_pCR))))
design <- model.matrix(~0 + grupos)
colnames(design) <- c("RD", "pCR")

fit <- lmFit(exprs(datos_TNBC), design)
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

plataforma <- getGEO(annotation(datos_TNBC), AnnotGPL = TRUE)
anotacion <- Table(plataforma)
relacion_sonda_gen <- data.frame(ID = anotacion$ID, GeneSymbol = anotacion$`Gene symbol`)

tabla_exportar <- data.frame(
  ID = ids_sondas[highlight],
  GeneSymbol = relacion_sonda_gen$GeneSymbol[match(ids_sondas[highlight], relacion_sonda_gen$ID)],
  LogFC = logFC_values[highlight],
  p_values = p_values[highlight],
  fdr_values = fdr_values[highlight]
)

print(tabla_exportar)
