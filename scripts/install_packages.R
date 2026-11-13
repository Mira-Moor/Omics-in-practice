# Package installation for BRB-seq analysis
# Run only when required packages are missing

.libPaths(c("R_library", .libPaths()))

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = "R_library", Ncpus = 1)
}

# Install Bioconductor packages if needed
for (pkg in c("DESeq2", "apeglm")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, lib = "R_library",
                         ask = FALSE, update = FALSE, Ncpus = 1)
  }
}

# Install CRAN packages if needed
for (pkg in c("pheatmap", "ggrepel")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, lib = "R_library", Ncpus = 1)
  }
}
