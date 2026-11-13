.libPaths(c("R_library", .libPaths()))

library(DESeq2)
library(data.table)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")

counts_df <- fread("data/gene_counts.tsv", data.table = FALSE, check.names = FALSE)
metadata <- fread("metadata/sample_metadata.tsv", data.table = FALSE)

gene_ids <- counts_df[[1]]
counts <- as.matrix(counts_df[, -1, drop = FALSE])
rownames(counts) <- gene_ids

mode(counts) <- "numeric"

stopifnot(!anyDuplicated(rownames(counts)))
stopifnot(!anyDuplicated(metadata$sample))
stopifnot(all(counts >= 0))
stopifnot(all(counts == round(counts)))
stopifnot(setequal(colnames(counts), metadata$sample))

metadata <- metadata[match(colnames(counts), metadata$sample), ]
stopifnot(identical(colnames(counts), metadata$sample))
rownames(metadata) <- metadata$sample

dim(counts)
summary(colSums(counts))
table(metadata$condition)
head(metadata)

 lib_df <- data.frame(sample = factor(metadata$sample_number, levels = 1:24),
                     library_size = colSums(counts),
                     condition = metadata$condition)

p_lib <- ggplot(lib_df, aes(sample, library_size, fill = condition)) +
  geom_col() +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Total raw counts", title = "Library sizes") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/library_sizes.pdf", p_lib, width = 8, height = 5)

metadata$condition <- factor(metadata$condition,
  levels = c("0h", "6h", "1d", "2d", "3d", "4d", "7d", "8d"))

metadata$replicate <- factor(metadata$replicate)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition
)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

dim(dds)

dds <- estimateSizeFactors(dds)
summary(sizeFactors(dds))

vsd <- vst(dds, blind = TRUE)
dim(assay(vsd))

pca_df <- plotPCA(vsd, intgroup = c("condition", "replicate"),
                  returnData = TRUE)

percent_var <- round(100 * attr(pca_df, "percentVar"))

pca_df$sample_number <- metadata[pca_df$name, "sample_number"]

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, shape = replicate)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = sample_number), size = 3, show.legend = FALSE) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  labs(title = "PCA of BRB-seq samples") +
  theme_bw()

ggsave("figures/PCA_samples.pdf", p_pca, width = 9, height = 6)

sample_dist <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dist)

sample_labels <- paste0("S", sprintf("%02d", metadata$sample_number),
                        "_", metadata$condition,
                        "_rep", metadata$replicate)

sample_annotation <- data.frame(
  time = metadata$condition,
  replicate = metadata$replicate
)
rownames(sample_dist_matrix) <- sample_labels
colnames(sample_dist_matrix) <- sample_labels
rownames(sample_annotation) <- sample_labels

pheatmap(sample_dist_matrix,
annotation_row = sample_annotation,
  annotation_col = sample_annotation,
  clustering_distance_rows = as.dist(sample_dist_matrix),
  clustering_distance_cols = as.dist(sample_dist_matrix),
  color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
  main = "Sample-to-sample distances (VST)",
  filename = "figures/sample_distance_heatmap.pdf",
  width = 9,
  height = 8)

# Select UBTF and 19 highly variable genes
ubtf_id <- "ENSG00000108312.17"

stopifnot(ubtf_id %in% rownames(vsd))

gene_variance <- apply(assay(vsd), 1, var)

top_genes <- head(setdiff(names(sort(gene_variance, decreasing = TRUE)),
                          ubtf_id), 19)

heatmap_genes <- c(ubtf_id, top_genes)

heatmap_mat <- assay(vsd)[heatmap_genes, , drop = FALSE]

colnames(heatmap_mat) <- sample_labels

# Gene symbol labels for the selected-gene heatmap
gene_symbols <- c(
  "ENSG00000108312.17" = "UBTF",
  "ENSG00000103168.18" = "TAF1C",
  "ENSG00000105849.7" = "POLR1F",
  "ENSG00000112592.16" = "TBP"
)


pheatmap(heatmap_mat,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_col = sample_annotation,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "Gene expression heatmap including UBTF",
  fontsize_row = 9,
  fontsize_col = 7,
  angle_col = 90,
  filename = "figures/UBTF_expression_heatmap.pdf",
  width = 12,
  height = 7)

# Section 8: Selected gene set heatmap

vsd_mat <- assay(vsd)

selected_genes <- scan("metadata/selected_genes.txt", what = character(), quiet = TRUE)

present <- intersect(selected_genes, rownames(vsd_mat))
missing <- setdiff(selected_genes, rownames(vsd_mat))

writeLines(missing, "results/selected_genes_not_found.txt")

if (length(present) < 2) stop("Fewer than two selected genes were found")

print(present)
print(missing)

# Section 8: Prepare the selected-gene z-score heatmap

heat_mat <- vsd_mat[present, , drop = FALSE]
colnames(heat_mat) <- sample_labels

heat_z <- t(scale(t(heat_mat)))

heat_z <- heat_z[complete.cases(heat_z), , drop = FALSE]

ann_col <- sample_annotation
# Replace Ensembl IDs with gene symbols
stopifnot(all(rownames(heat_z) %in% names(gene_symbols)))

rownames(heat_z) <- unname(gene_symbols[rownames(heat_z)])

# Section 8: Selected-gene row z-score heatmap

pdf("figures/selected_gene_heatmap_zscore.pdf", width = 8, height = 7)

pheatmap(heat_z,
  annotation_col = ann_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = TRUE,
  main = "Selected genes: row z-scores")

dev.off()

# Section 9: Heatmap in experimental order

sample_order <- sample_labels[order(metadata$condition,
                                    metadata$replicate)]

heat_z_ordered <- heat_z[, sample_order, drop = FALSE]

pdf("figures/selected_gene_heatmap_ordered.pdf",
    width = 8, height = 7)

pheatmap(heat_z_ordered,
  annotation_col = ann_col[sample_order, , drop = FALSE],
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  main = "Selected genes by experimental group")

dev.off()

# Section 9: Unscaled VST expression heatmap

heat_mat_vst <- heat_mat

stopifnot(all(rownames(heat_mat_vst) %in% names(gene_symbols)))

rownames(heat_mat_vst) <- unname(gene_symbols[rownames(heat_mat_vst)])

stopifnot(identical(colnames(heat_mat_vst), rownames(ann_col)))

pdf("figures/selected_gene_heatmap_VST.pdf",
    width = 8, height = 7)

pheatmap(heat_mat_vst,
  annotation_col = ann_col,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = TRUE,
  main = "Selected genes: VST expression")

dev.off()

# Section 10: Differential expression analysis with DESeq2

dds <- DESeq(dds)

resultsNames(dds)

# Section 10: Differential expression at 6h versus 0h

res <- results(
  dds,
  contrast = c("condition", "6h", "0h"),
  alpha = 0.05
)

summary(res)

# Examine UBTF specifically
print(res[ubtf_id, ])

# Section 10: Shrink log2 fold changes (6h vs 0h)

res_shrunk <- lfcShrink(
  dds,
  coef = "condition_6h_vs_0h",
  type = "apeglm"
)

# Examine the shrunken UBTF result
print(res_shrunk[ubtf_id, ])

# Section 10: Export differential expression results

res_df <- as.data.frame(res_shrunk)

res_df$gene_id <- rownames(res_df)

res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

fwrite(res_df,
       "results/DE_6h_vs_0h.tsv",
       sep = "\t",
       na = "NA")

# Section 11: Classify differential expression results

padj_cutoff <- 0.05
lfc_cutoff <- 1

res_df$status <- "Not significant"

res_df$status[
  !is.na(res_df$padj) &
  res_df$padj < padj_cutoff &
  res_df$log2FoldChange >= lfc_cutoff
] <- "Higher at 6h"

res_df$status[
  !is.na(res_df$padj) &
  res_df$padj < padj_cutoff &
  res_df$log2FoldChange <= -lfc_cutoff
] <- "Lower at 6h"

print(table(res_df$status))

# Section 11: MA plot for 6h versus 0h

pdf("figures/MA_6h_vs_0h.pdf", width = 7, height = 5)

plotMA(res_shrunk,
       alpha = 0.05,
       ylim = c(-5, 5),
       main = "6h versus 0h")

dev.off()

# Section 12: Volcano plot for 6h versus 0h

plot_df <- res_df

plot_df$minus_log10_padj <- -log10(
  pmax(plot_df$padj, .Machine$double.xmin)
)

plot_df$minus_log10_padj[is.na(plot_df$padj)] <- NA_real_

label_pool <- subset(plot_df,
  status != "Not significant" & !is.na(padj))

label_genes <- head(label_pool[order(label_pool$padj), ], 12)

p_volcano <- ggplot(plot_df,
  aes(log2FoldChange, minus_log10_padj, color = status)) +
  geom_point(alpha = 0.65, size = 1.4, na.rm = TRUE) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
             linetype = 2) +
  geom_hline(yintercept = -log10(padj_cutoff),
             linetype = 2) +
  scale_color_manual(values = c(
    "Lower at 6h" = "#377EB8",
    "Not significant" = "grey70",
    "Higher at 6h" = "#E41A1C")) +
  geom_text_repel(data = label_genes,
                  aes(label = gene_id),
                  size = 3, max.overlaps = Inf,
                  show.legend = FALSE) +
  labs(title = "6h versus 0h",
       x = "Shrunken log2 fold change",
       y = expression(-log[10](adjusted~p)),
       color = NULL) +
  theme_bw()

ggsave("figures/volcano_6h_vs_0h.pdf",
       p_volcano, width = 7.5, height = 6)

ggsave("figures/volcano_6h_vs_0h.png",
       p_volcano, width = 7.5, height = 6, dpi = 300)

# Section 13: Validate UBTF expression across time points

count_df <- plotCounts(
  dds,
  gene = ubtf_id,
  intgroup = "condition",
  returnData = TRUE
)

p_count <- ggplot(count_df, aes(x = condition, y = count)) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 3
  ) +
  scale_y_log10() +
  labs(
    title = "UBTF expression across time points",
    x = "Time point",
    y = "Normalized count (log scale)"
  ) +
  theme_bw()

ggsave(
  "figures/counts_UBTF.pdf",
  p_count,
  width = 7,
  height = 5
)

# Section 14: Final reproducibility exports

# Save PCA coordinates and sample information
fwrite(
  pca_df,
  "results/PCA_coordinates.tsv",
  sep = "\t",
  na = "NA"
)

# Update the DESeq2 table to include status classifications
fwrite(
  res_df,
  "results/DE_6h_vs_0h.tsv",
  sep = "\t",
  na = "NA"
)
