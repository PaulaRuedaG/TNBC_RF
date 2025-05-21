library(Biobase)
library(limma)    

datos <- read.csv("C:/Users/pasan/Documents/Maestria/Trabajo_de_grado/Analisis_conjunto_DEGs/fRMA_normalized_expression.csv", row.names = 1)  # Usa row.names = 1 si los genes están en la primera columna

# Limpiar los nombres de las columnas
colnames(datos) <- gsub("_.*", "", colnames(datos)) 

#Cargar nombres de las muestras
muestras_TNBC_RD1 <- readLines("C:/Users/pasan/Documents/Maestria/Trabajo_de_grado/Analisis_conjunto_DEGs/muestras_TNBC_RD_GSE25066.txt") 
muestras_TNBC_pCR1 <- readLines("C:/Users/pasan/Documents/Maestria/Trabajo_de_grado/Analisis_conjunto_DEGs/muestras_TNBC_pCR_GSE25066.txt") 
muestras_TNBC_RD2 <- readLines("C:/Users/pasan/Documents/Maestria/Trabajo_de_grado/Analisis_conjunto_DEGs/muestras_TNBC_RD_GSE20194.txt") 
muestras_TNBC_pCR2 <- readLines("C:/Users/pasan/Documents/Maestria/Trabajo_de_grado/Analisis_conjunto_DEGs/muestras_TNBC_pCR_GSE20194.txt") 
muestras_TNBC_RD <- c(muestras_TNBC_RD1,muestras_TNBC_RD2)
muestras_TNBC_pCR <- c(muestras_TNBC_pCR1,muestras_TNBC_pCR2)
muestras_seleccionadas <- c(muestras_TNBC_RD,muestras_TNBC_pCR)


# Crear un vector de grupos basado en los nombres de las columnas
grupo <- ifelse(colnames(datos) %in% muestras_TNBC_RD, "RD", 
                ifelse(colnames(datos) %in% muestras_TNBC_pCR, "pCR", NA))  # NA si alguna no coincide

# Crear el dataframe de metadata
metadata <- data.frame(Sample = colnames(datos), Group = grupo)
metadata <- metadata[!is.na(metadata$Group), ]  # Eliminar filas con NA en la columna "Group"

# Asegurar que metadata solo contenga las muestras que están en data
metadata <- metadata[metadata$Sample %in% colnames(datos), ]

# Asegurar que data solo contenga las muestras que están en metadata
datos <- datos[, colnames(datos) %in% metadata$Sample]

# Ordenar metadata para que coincida con el orden de las columnas en data
metadata <- metadata[match(colnames(datos), metadata$Sample), ]

# Verificar que los nombres de las muestras coincidan
stopifnot(all(colnames(datos) == metadata$Sample))

# Convertir metadata en AnnotatedDataFrame
rownames(metadata) <- metadata$Sample  # Usar Sample como rownames
metadata$Sample <- NULL  # Eliminar la columna redundante
metadata <- new("AnnotatedDataFrame", data = metadata)  # Creación correcta

# Convertir la matriz de expresión en AssayData
assay_data <- assayDataNew("environment", exprs = as.matrix(datos))

# Crear el objeto ExpressionSet
expSet <- ExpressionSet(assayData = as.matrix(datos), phenoData = metadata)

# Verificar la estructura del objeto ExpressionSet
expSet
pData(expSet)  # Ver metadata dentro del ExpressionSet
exprs(expSet)  # Ver datos de expresión

realizar_bootstrap_DEGs <- function(muestras_TNBC_RD, muestras_TNBC_pCR, expSet, num_iter = num_iteraciones, logFC_threshold, p_value_threshold) {
  
  # Crear un data frame vacío para almacenar los resultados de cada iteración
  lista_acumulada_DEGs <- list()
  
  for (i in 1:num_iter) {
    cat("Iteración: ", i, "\n")
    
    # Tomar muestras aleatorias con reemplazo para RD y pCR
    muestras_RD_bootstrap <- sample(muestras_TNBC_RD2, 50, replace = TRUE)
    muestras_pCR_bootstrap <- sample(muestras_TNBC_pCR2, 50, replace = TRUE)
    muestras_bootstrap <- c(muestras_RD_bootstrap, muestras_pCR_bootstrap) 
    datos_bootstrap <- expSet[, muestras_bootstrap]
    
    # Definir el factor para los grupos RD y pCR
    grupos_bootstrap <- factor(c(rep("RD", length(muestras_RD_bootstrap)), rep("pCR", length(muestras_pCR_bootstrap))))
    
    # Crear la matriz de diseño
    design_bootstrap <- model.matrix(~0 + grupos_bootstrap)
    colnames(design_bootstrap) <- c("RD", "pCR")
    
    # Ajuste del modelo lineal usando limma
    fit_bootstrap <- lmFit(exprs(datos_bootstrap), design_bootstrap)
    contraste_bootstrap <- makeContrasts(RD_vs_pCR = RD - pCR, levels = design_bootstrap)
    fit_bootstrap2 <- contrasts.fit(fit_bootstrap, contrasts = contraste_bootstrap)
    
    # Ajuste estadístico
    fit_bootstrap2 <- eBayes(fit_bootstrap2)
    
    # Extraer la tabla de resultados completa con ajuste por Benjamini-Hochberg (FDR)
    tT_bootstrap <- topTable(fit_bootstrap2, coef = 1, adjust = "fdr", number = Inf)
    
    # Filtrar los genes que cumplen con el umbral de logFC y p-value ajustado
    DEGs_filtrados <- tT_bootstrap[abs(tT_bootstrap$logFC) > logFC_threshold & tT_bootstrap$adj.P.Val < p_value_threshold, ]
    
    # Almacenar solo los nombres de los genes de esta iteración
    lista_acumulada_DEGs[[i]] <- rownames(DEGs_filtrados)
  }
  
  # Combinar todos los DEGs de todas las iteraciones en un solo vector
  DEGs_combinados <- unlist(lista_acumulada_DEGs)
  #print(DEGs_combinados)
  
  # Retornar la lista de DEGs acumulados
  return(DEGs_combinados)
}

#---------------------------------------------------
# Ejecutar el bootstrap
#---------------------------------------------------

num_iteraciones <- 999  # Número de iteraciones del bootstrap
logFC_threshold <- 1  # Umbral para logFC
p_value_threshold <- 0.05  # Umbral de significancia estadística

# Ejecutar la función y almacenar todos los resultados
resultados_bootstrap <- realizar_bootstrap_DEGs(muestras_TNBC_RD, muestras_TNBC_pCR, expSet, num_iteraciones, logFC_threshold, p_value_threshold)

#---------------------------------------------------
# Almacenar los DEGs
#---------------------------------------------------

DEGs_combinados <- unlist(resultados_bootstrap)
tabla_frecuencias_DEGs <- as.data.frame(table(DEGs_combinados))
tabla_frecuencias_DEGs <- tabla_frecuencias_DEGs[order(-tabla_frecuencias_DEGs$Freq), ]
head(tabla_frecuencias_DEGs)

# Seleccionar los DEGs presentes en n% de las iteraciones
umbral_frecuencia <- 0.2  # Umbral de frecuencia (30%)
genes_destacados <- tabla_frecuencias_DEGs$DEGs_combinados[tabla_frecuencias_DEGs$Freq >= (num_iteraciones * umbral_frecuencia)]

#---------------------------------------------------
# Graficar el volcano plot 
#---------------------------------------------------

#Generar un volcano plot completo RD vs pCR

# Definir grupos y diseño experimental
grupos <- factor(c(rep("RD", length(muestras_TNBC_RD)), rep("pCR", length(muestras_TNBC_pCR))))
design <- model.matrix(~0 + grupos)
colnames(design) <- c("RD", "pCR")

# Ajuste del modelo lineal con limma
fit <- lmFit(exprs(expSet), design)
contraste <- makeContrasts(RD_vs_pCR = RD - pCR, levels = design)
fit2 <- contrasts.fit(fit, contrasts = contraste)
fit2 <- eBayes(fit2)

# Extraer resultados con FDR ajustado
resultados <- topTable(fit2, coef=1, number=Inf, adjust="fdr")

# Usar directamente los IDs de las sondas como etiquetas
ids_sondas <- rownames(resultados)

# Selección de genes destacados usando IDs
genes_destacados_character <- as.character(genes_destacados)
highlight <- which(ids_sondas %in% genes_destacados_character)

# Extraer valores para el volcano plot
logFC_values <- resultados$logFC
p_values <- resultados$P.Value
fdr_values <- resultados$adj.P.Val

# Generar el volcano plot usando IDs
volcanoplot(fit2, coef=1, main="Volcano Plot RD vs pCR", pch=20,
            highlight=highlight, col=ifelse(1:nrow(fit2) %in% highlight, "red", "black"),
            xlab = "Log2 Fold Change", ylab = "-Log10(Adjusted P-value)")

# Agregar los IDs de las sondas destacadas en el gráfico
highlight_logFC <- logFC_values[highlight]
highlight_p_value <- -log10(p_values[highlight])

# Agregar las etiquetas de los IDs de las sondas
text(highlight_logFC, highlight_p_value, labels = ids_sondas[highlight], 
     pos = 3, col = "red", cex = 0.7)
legend("topright", legend = c("Probes destacados"), col = "red", pch = 20)

# Imprimir los IDs de las sondas destacadas al final
cat("IDs de las sondas resaltadas en el volcano plot:\n")
candidatos <- ids_sondas[highlight]
print(candidatos)

# Obtener la tabla de anotación de GEO
relacion_sonda_gen <- read.csv("C:/Users/pasan/Documents/Maestria/Trabajo_de_grado/Analisis_conjunto_DEGs/relacion_sonda_gen.csv")

# Unir los datos con la tabla de resultados usando el ID de la sonda
tabla_exportar <- data.frame(
  ID = ids_sondas,
  GeneSymbol = relacion_sonda_gen$GeneSymbol[match(ids_sondas, relacion_sonda_gen$ID)],
  LogFC = logFC_values,
  p_values = p_values,
  fdr_values = fdr_values
)

print(tabla_exportar)

write.csv(tabla_exportar, file = "genes_destacados_volcano4.csv", row.names = FALSE)
cat("Tabla exportada exitosamente a 'genes_destacados_volcano.csv'\n")

#--------------------------------------------------------------------------------------
#genes_interes <- c("38241_at", "216685_s_at", "213033_s_at", "201656_at", "211363_s_at", "219850_s_at", "205833_s_at", "209686_at", "204775_at", "1438_at", "204351_at", "215551_at", "209290_s_at", "206686_at", "208712_at","201755_at", "205225_at", "204107_at", "220624_s_at", "205240_at", "217163_at", "219051_x_at", "211627_x_at", "222031_at", "217190_x_at", "211234_x_at", "201438_at", "201504_s_at", "208711_s_at", "202310_s_at", "211466_at", "219197_s_at", "210683_at", "211233_x_at", "202273_at", "219970_at", "204886_at", "205478_at", "215552_s_at", "211235_s_at", "219654_at", "209644_x_at")

#genes_interes <- c("204107_at","205833_s_at","213033_s_at","206686_at","204886_at","222031_at","219051_x_at","208711_s_at","208712_at","201504_s_at","209686_at","201656_at","219850_s_at","211363_s_at","219654_at","211466_at","201438_at","205478_at","201755_at","216685_s_at","204775_at","205240_at","202273_at","202310_s_at","220624_s_at","210683_at","204351_at","219970_at","1438_at","209290_s_at","209644_x_at","38241_at","205225_at", "211233_x_at", "211234_x_at", "211235_s_at", "211627_x_at", "215551_at", "215552_s_at", "217163_at", "217190_x_at","219197_s_at","212595_s_at","213826_s_at","202648_at","213350_at","221943_x_at")

# def genes_interes <- c("204107_at","205833_s_at","213033_s_at","206686_at","204886_at","222031_at","219051_x_at","208711_s_at","208712_at","201504_s_at","209686_at","201656_at","219850_s_at","211363_s_at","219654_at","211466_at","201438_at","205478_at","201755_at","216685_s_at","204775_at","205240_at","202273_at","202310_s_at","220624_s_at","210683_at","204351_at","219970_at","1438_at","209290_s_at","209644_x_at","38241_at","205225_at","211233_x_at","211234_x_at","211235_s_at","211627_x_at","215551_at","215552_s_at","217163_at","217190_x_at","219197_s_at")


mtrx_expresion=exprs(expSet)
mtrx_expresion_filtrada <- as.data.frame(t(mtrx_expresion[, muestras_seleccionadas, drop = FALSE]))
mtrx_expresion_filtrada$grupos <- grupos
degs_filtrados<- mtrx_expresion_filtrada[, c(genes_interes, "grupos"), drop = FALSE]
write.csv(degs_filtrados, "datos_ML_Conjuntos_7.csv", row.names = TRUE)

