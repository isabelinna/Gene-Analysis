library(GEOquery)
library(limma)

# 1 - Leer el conjunto de datos
gset <- getGEO("GSE17536", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# 2 - Seleccionar solo las sondas que cuentan con symbolo de gen
id=gset@featureData@data[["ID"]]
symbol=data.frame(gset@featureData@data[["Gene symbol"]])
row.names(symbol) = id
gene_id = id[symbol != ""]
gene_symbol = data.frame(symbol[symbol != ""])
row.names(gene_symbol) = gene_id

# 3 - Obtiene los valores de expresión
ex <- exprs(gset)

# 4 - Seleccionar 10 muestras de pacientes en etapa 1 y 10 en etapa 3 
ex <- ex[,c(2,4,7,15,23,41,46,49,58,79,26,35,50,51,55,65,92,95,104,110)]
ex <- ex[row.names(gene_symbol),]

# 5 - Normalizacion
raw_means = apply(ex,2,mean,trim=0.02)

microarray_norm = sweep(ex, 2, raw_means, "/") * 100

# 6 - Medias
cancer_estado1_mean = rowMeans(microarray_norm[,1:10])
cancer_estado3_mean = rowMeans(microarray_norm[,11:20])

microarray_means = data.frame(cancer_estado1_mean, cancer_estado3_mean)

# 7 - Proporciones
cancer_ratios = microarray_means$cancer_estado3_mean / microarray_means$cancer_estado1_mean

microarray_ratios = data.frame(cancer_ratios)

row.names(microarray_ratios) = row.names(ex)


# 8 - Cambio a log2
microarray_norm_log2 = log2(microarray_norm)
microarray_means_log2 = log2(microarray_means)
microarray_ratios_log2 = log2(microarray_ratios)

microarray_ratios_log2_v = unlist(microarray_ratios_log2)
names(microarray_ratios_log2_v) = row.names(microarray_ratios)

# 9 - t-test
get_pvalue <- function(values, idx1, idx2) {
  return(t.test(values[idx1], values[idx2])$p.value)
}

cancer_p = apply(microarray_norm, 1, get_pvalue, 1:10, 11:20)


# 10 - Filtrar genes con valor de significancia p < 0.05
filtered_cancer_p = cancer_p[cancer_p < 0.05]

# 11 - Filtrar genes con tamano de efecto grande log2 > 0.38 y log2 < -0.38
filtered_cancer_log2 = c(microarray_ratios_log2_v[microarray_ratios_log2_v > 0.38], 
                         microarray_ratios_log2_v[microarray_ratios_log2_v < -0.38])

# 12 - Proporciones en log2 de sondas con p<0.05
sondas_log2_p = microarray_ratios_log2_v[names(filtered_cancer_p)]
slp = data.frame(sondas_log2_p)

# 13 - Top 10 genes que se sobreexpresan y subexpresan
sobreexp = sort(sondas_log2_p, decreasing = TRUE)[1:10]
subexp = sort(sondas_log2_p)[1:10]

genes_sobre_sub = data.frame("Genes sobre" = gene_symbol[names(sobreexp),], 
                             "valor exp sobre" = sobreexp, "Genes sub" = gene_symbol[names(subexp),], 
                             "valor exp sub" = subexp) 
row.names(genes_sobre_sub) = c(1:10)

# 14 - Genes con valor de significancia p < 0.5 y tamaño de efecto grande.
genes_final_names = intersect(names(filtered_cancer_p), names(filtered_cancer_log2))
genes_final = data.frame("Gen" = gene_symbol[genes_final_names,], 
                         "Expresion en Log2" = microarray_ratios_log2_v[genes_final_names])

# 15 - Generación de gráficas.
# Mapas de calor y dendogramas
microarray_selection = microarray_means[genes_final_names,]

medias = rowMeans(microarray_selection)
devs = apply(microarray_selection, 1, sd)

centered_microarray_selection = sweep(microarray_selection, 1, medias)
centered_microarray_selection = sweep(centered_microarray_selection, 1, devs, "/")

names(centered_microarray_selection) = c("cancer1", "cancer3")
hclustering = hclust(dist(centered_microarray_selection))

plot(hclustering)

names(microarray_selection) = c("E1", "E2")
heatmap(as.matrix(microarray_selection), Colv = NA, 
        main = "Mapa de calor: Valores de expresión")

# Gráfica de dispersión
plot(microarray_means_log2$cancer_estado1_mean,
     microarray_means_log2$cancer_estado3_mean,
     xlim = c(0, 9), ylim = c(0, 9),
     xaxt = "n", yaxt = "n",
     main = "Gráfica de dispersión: Expresión en Etapa 3 vs Etapa 1",
     ylab = "Etapa3 (expresión en log2)",
     xlab = "Etapa1 (expresión en log2)")

axis(1, at = seq(0:9))
axis(2, at = seq(0:9))
abline(lm(microarray_means_log2$cancer_estado1_mean ~ microarray_means_log2$cancer_estado3_mean), 
       col="red")

# Gráfica R-I
plot(microarray_means_log2$cancer_estado1_mean + microarray_means_log2$cancer_estado3_mean,
     microarray_means_log2$cancer_estado1_mean - microarray_means_log2$cancer_estado3_mean,
     main = "Gráfica R-I: Expresión en Etapa 3 vs Etapa 1",
     xlab ="log2(Etapa 3 * Etapa 1)",
     ylab = "log2(Etapa 3 / Etapa 1)")

# Gráfica volcán
colores = rep(1, length(cancer_p))
colores[cancer_p < 0.05 & microarray_ratios_log2 < -0.38] = 2
colores[cancer_p < 0.05 & microarray_ratios_log2 > 0.38] = 3

plot(microarray_ratios_log2_v, cancer_p, col=colores, log = "y", ylim = rev(range(cancer_p)),
     main = "Gráfica de volcán: Sobre expresión y sub expresión",
     xlab = "Proporción de expresion en log2: Etapa 3 vs Etapa 1",
     ylab = "p-value")


