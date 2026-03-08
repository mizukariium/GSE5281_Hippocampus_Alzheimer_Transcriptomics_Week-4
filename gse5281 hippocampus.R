BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("enrichplot")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("dplyr")

library(BiocManager)
library(GEOquery)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(pheatmap)
library(dplyr)

library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(hgu133a.db)
library(AnnotationDbi)
library(umap)

gset <- getGEO("GSE5281", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
View(gset)
library(GEOquery)
gset <- getGEO("GSE5281", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
View(gset)
gset
pData(gset)
ex <- exprs(gset)
View(ex)

library(GEOquery)
gset <- getGEO("GSE5281", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
pheno <- pData(gset)
head(pheno)
colnames(pheno)
"brain region:"
"disease state:"
hippo_idx <- grep("hippocampus", pheno$`brain region:`, ignore.case = TRUE)
gset_hippo <- gset[, hippo_idx]
pheno_hippo <- pData(gset_hippo)
table(pheno_hippo$`brain region:`)
colnames(pheno)
head(pheno$characteristics_ch1)
head(pheno$characteristics_ch1.1)
head(pheno$characteristics_ch1.2)
head(pheno$title)
hippo_idx <- grep("^HIP", pheno$title)
length(hippo_idx)
gset_hippo <- gset[, hippo_idx]
pheno_hippo <- pData(gset_hippo)
table(pheno_hippo$title)
group <- ifelse(
  grepl("control", pheno_hippo$title, ignore.case = TRUE),
  "Control",
  "AD"
)
group <- factor(group, levels = c("Control", "AD"))
table(group)
expr <- exprs(gset_hippo)
dim(expr)

View(ex)
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}
group_info <- pData(gset)[["source_name_ch1"]]
groups <- make.names(group_info)
gset$group <- factor(groups)

library(GEOquery)
library(limma)

gset <- getGEO("GSE5281", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
ex <- exprs(gset)

qx <- as.numeric(quantile(ex,
                          c(0, 0.25, 0.5, 0.75, 0.99, 1),
                          na.rm = TRUE))

LogTransform <- (qx[5] > 100) ||
  (qx[6] - qx[1] > 50 && qx[2] > 0)

if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

pheno <- pData(gset)

hippo_idx <- grep(
  "hippocampus|HIP",
  apply(pheno, 1, paste, collapse = " "),
  ignore.case = TRUE
)

gset <- gset[, hippo_idx]
ex <- ex[, hippo_idx]

pheno <- pData(gset)

dim(ex)

group_info <- apply(pheno, 1, paste, collapse = " ")

groups <- ifelse(
  grepl("control", group_info, ignore.case = TRUE),
  "Control",
  "AD"
)

gset$group <- factor(groups)

nama_grup <- levels(gset$group)
print(nama_grup)


design <- model.matrix(~0 + gset$group)
colnames(design) <- levels(gset$group)
design

grup_ad <- "AD"
grup_control <- "Control"

contrast_formula <- paste(grup_ad, "-", grup_control)

print(paste("Kontras yang dianalisis:",
            contrast_formula))

contrast.matrix <- makeContrasts(
  ADvsControl = AD - Control,
  levels = design
)

fit <- lmFit(ex, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

deg <- topTable(fit2,
                adjust.method = "BH",
                number = Inf)

head(deg)

fit <- lmFit(ex, design)

contrast_matrix <- makeContrasts(
  contrasts = contrast_formula,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)

fit2 <- eBayes(fit2)

topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.01
)

head(topTableResults)

View(topTableResults)

probe_ids <- rownames(topTableResults)

gene_annotation <- AnnotationDbi::select(
  hgu133a.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)
topTableResults$PROBEID <- rownames(topTableResults)

topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)
View(topTableResults)

head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])

#1 BOXPLOT DISTRIBUSI EKSPRESI
group_colors <- as.numeric(gset$group)

boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot Distribusi Nilai Ekspresi Hippocampus",
  ylab = "Expression Value (log2)"
)

legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.8
)

#2 DENSITY PLOT DISTRIBUSI EKSPRESI
library(ggplot2)

expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen Hippocampus",
    x = "Expression Value (log2)",
    y = "Density"
  )

#3 UMAP PLOT (VISUALISASI CLUSTERING SAMPLE)
library(umap)

umap_input <- t(ex)

umap_result <- umap(umap_input)

umap_df <- data.frame(
  UMAP1 = umap_result$layout[,1],
  UMAP2 = umap_result$layout[,2],
  Group = gset$group
)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot Sampel Hippocampus (AD vs Control)",
    x = "UMAP 1",
    y = "UMAP 2"
  )

#4 VOLCANO PLOT
volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)

volcano_data$status <- "NO"

volcano_data$status[
  volcano_data$logFC > 1 &
    volcano_data$adj.P.Val < 0.05
] <- "UP"

volcano_data$status[
  volcano_data$logFC < -1 &
    volcano_data$adj.P.Val < 0.05
] <- "DOWN"

ggplot(volcano_data,
       aes(x = logFC,
           y = -log10(adj.P.Val),
           color = status)) +
  
  geom_point(alpha = 0.6) +
  
  scale_color_manual(values = c(
    "DOWN" = "blue",
    "NO" = "grey",
    "UP" = "red"
  )) +
  
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  
  theme_minimal() +
  
  ggtitle("Volcano Plot Differentially Expressed Genes Hippocampus AD vs Control")


#5 HEATMAP TOP 50 DEG
library(pheatmap)

topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]

top50 <- head(topTableResults, 50)

mat_heatmap <- ex[top50$PROBEID, ]

gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID,
  top50$SYMBOL
)

rownames(mat_heatmap) <- gene_label

# hapus NA
mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]

# hapus varians nol
gene_variance <- apply(mat_heatmap, 1, var)

mat_heatmap <- mat_heatmap[
  gene_variance > 0,
]

annotation_col <- data.frame(
  Group = gset$group
)

rownames(annotation_col) <- colnames(mat_heatmap)

pheatmap(
  mat_heatmap,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes Hippocampus AD vs Control"
)

#6 SIMPAN HASIL DEG KE FILE CSV
write.csv(
  topTableResults,
  "Hasil_GSE5281_Hippocampus_DEG.csv"
)

message("Analisis GSE5281 Hippocampus selesai. File hasil disimpan.")

#7 GO ENRICHMENT ANALYSIS
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
##7.1
deg_sig <- topTableResults[
  topTableResults$adj.P.Val < 0.05 &
    !is.na(topTableResults$SYMBOL),
]
gene_symbols <- deg_sig$SYMBOL
##7.2
gene_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_symbols,
  columns = "ENTREZID",
  keytype = "SYMBOL"
)

gene_entrez <- gene_entrez$ENTREZID

gene_entrez <- na.omit(gene_entrez)
##7.3
ego <- enrichGO(
  gene = gene_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)
##7.4
dotplot(
  ego,
  showCategory = 20,
  title = "GO Biological Process Enrichment - Hippocampus AD vs Control"
)

#8 KEGG PATHWAY ENRICHMENT
ekegg <- enrichKEGG(
  gene = gene_entrez,
  organism = "hsa",
  pvalueCutoff = 0.05
)
dotplot(
  ekegg,
  showCategory = 20,
  title = "KEGG Pathway Enrichment - Hippocampus AD vs Control"
)


write.csv(
  as.data.frame(ego),
  "GO_Enrichment_GSE5281.csv"
)

write.csv(
  as.data.frame(ekegg),
  "KEGG_Enrichment_GSE5281.csv"
)

expr_matrix <- exprs(gset)
write.csv(expr_matrix, "GSE5281_hippocampus_expression.csv")
pheno_data <- pData(gset)
write.csv(pheno_data, "GSE5281_metadata.csv")

