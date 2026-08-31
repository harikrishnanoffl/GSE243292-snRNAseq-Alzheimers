# ================================================================
# GSE243292 Alzheimer's Disease snRNA-seq Analysis
# Complete Seurat 5 Workflow
# ================================================================
#
# Python:
#   Used ONLY to convert the original H5AD file into sparse
#   Matrix Market files.
#
# R / Seurat:
#   Performs the complete biological analysis:
#   1. Load converted data
#   2. Quality control
#   3. Normalization
#   4. Highly variable genes
#   5. Scaling
#   6. PCA
#   7. Clustering
#   8. UMAP
#   9. t-SNE
#   10. Cell-type annotation
#   11. Marker gene identification
#   12. Differential expression
#   13. GO enrichment
#

# ================================================================
# SECTION 1 — SET PROJECT ROOT
# ================================================================
#
# IMPORTANT:
# The code uses relative paths after setting the working directory.
#
# For GitHub, another user can simply open the repository and set
# the working directory to the repository root.
#
# ================================================================

setwd("D:/nsg consultaion")

getwd()


# ================================================================
# SECTION 2 — INSTALL REQUIRED PACKAGES
# ================================================================
#
# Run this section only once on a new computer.
#
# ================================================================

cran_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "ggplot2",
  "patchwork",
  "dplyr",
  "data.table",
  "R.utils",
  "jsonlite",
  "tidyr"
)

installed_packages <- rownames(
  installed.packages()
)

for (pkg in cran_packages) {
  
  if (!pkg %in% installed_packages) {
    
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}


# Load packages

suppressPackageStartupMessages({
  
  library(Seurat)
  
  library(SeuratObject)
  
  library(Matrix)
  
  library(ggplot2)
  
  library(patchwork)
  
  library(dplyr)
  
  library(data.table)
  
  library(R.utils)
  
  library(jsonlite)
  library(tidyr)
  
})


# Check versions

R.version.string

packageVersion("Seurat")

packageVersion("SeuratObject")


# ================================================================
# SECTION 3 — CREATE PROJECT DIRECTORIES
# ================================================================

dir.create(
  "results",
  showWarnings = FALSE
)

dir.create(
  "figures",
  showWarnings = FALSE
)

dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)


# ================================================================
# SECTION 4 — DEFINE INPUT FILES
# ================================================================
#
# These files were created by:
#
# python/01_convert_h5ad.py
#
# ================================================================

processed_dir <- "data/processed"

matrix_file <- file.path(
  processed_dir,
  "matrix.mtx.gz"
)

features_file <- file.path(
  processed_dir,
  "features.tsv.gz"
)

barcodes_file <- file.path(
  processed_dir,
  "barcodes.tsv.gz"
)

cell_metadata_file <- file.path(
  processed_dir,
  "cell_metadata.tsv"
)

gene_metadata_file <- file.path(
  processed_dir,
  "gene_metadata.tsv"
)


# Check files

required_files <- c(
  matrix_file,
  features_file,
  barcodes_file,
  cell_metadata_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  stop(
    paste(
      "The following files are missing:",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}


# ================================================================
# SECTION 5 — LOAD SPARSE EXPRESSION MATRIX
# ================================================================
#
# Python converted the H5AD matrix into:
#
# matrix.mtx.gz
# features.tsv.gz
# barcodes.tsv.gz
#
# ReadMtx() loads the matrix without converting it to a dense
# matrix, which is important for a 12 GB RAM laptop.
#
# ================================================================

counts <- ReadMtx(
  
  mtx = matrix_file,
  
  cells = barcodes_file,
  
  features = features_file,
  
  cell.column = 1,
  
  feature.column = 2,
  
  cell.sep = "\t",
  
  feature.sep = "\t",
  
  unique.features = TRUE,
  
  strip.suffix = FALSE
  
)


# Confirm sparse matrix

if (!inherits(
  counts,
  "sparseMatrix"
)) {
  
  stop(
    "Expression matrix is not sparse."
  )
}


cat(
  "Number of genes:",
  nrow(counts),
  "\n"
)

cat(
  "Number of cells:",
  ncol(counts),
  "\n"
)


# ================================================================
# SECTION 6 — LOAD CELL METADATA
# ================================================================

cell_metadata <- read.delim(
  
  cell_metadata_file,
  
  header = TRUE,
  
  sep = "\t",
  
  check.names = FALSE,
  
  stringsAsFactors = FALSE
  
)


# Check cell_id

if (!"cell_id" %in%
    colnames(cell_metadata)) {
  
  stop(
    "cell_metadata.tsv does not contain cell_id."
  )
}


# Set cell IDs as row names

cell_metadata$cell_id <- as.character(
  cell_metadata$cell_id
)

rownames(cell_metadata) <- (
  cell_metadata$cell_id
)


# ================================================================
# SECTION 7 — VERIFY CELL IDs
# ================================================================

matrix_cells <- colnames(counts)

metadata_cells <- rownames(
  cell_metadata
)


# Cells missing from metadata

missing_metadata <- setdiff(
  
  matrix_cells,
  
  metadata_cells
  
)


if (length(missing_metadata) > 0) {
  
  stop(
    paste(
      "Some cells do not have metadata.",
      "Number missing:",
      length(missing_metadata)
    )
  )
}


# Reorder metadata to match expression matrix

cell_metadata <- cell_metadata[
  matrix_cells,
  ,
  drop = FALSE
]


# Final verification

if (!identical(
  
  rownames(cell_metadata),
  
  colnames(counts)
  
)) {
  
  stop(
    "Cell metadata and expression matrix do not match."
  )
}


cat(
  "Cell IDs successfully matched.\n"
)


# ================================================================
# SECTION 8 — CREATE SEURAT OBJECT
# ================================================================
#
# CreateSeuratObject() creates the main Seurat object.
#
# min.cells = 10:
#   Keep genes detected in at least 10 cells.
#
# min.features = 200:
#   Remove extremely low-quality cells.
#
# More detailed QC is performed later.
#
# ================================================================

seurat_obj <- CreateSeuratObject(
  
  counts = counts,
  
  assay = "RNA",
  
  project = "GSE243292",
  
  meta.data = cell_metadata,
  
  min.cells = 10,
  
  min.features = 200
  
)


# Remove standalone matrix from memory

rm(counts)

rm(cell_metadata)

gc()


# Inspect object

seurat_obj


# ================================================================
# SECTION 9 — MAP THE 15 GSE243292 SAMPLES
# ================================================================
#
# Sample groups:
#
# Normal:
#   GSM7782916
#   GSM7782917
#
# Pathological Aging:
#   GSM7782918
#   GSM7782919
#   GSM7782920
#   GSM7782921
#   GSM7782922
#
# Alzheimer's Disease:
#   GSM7782923
#   GSM7782924
#   GSM7782925
#   GSM7782926
#   GSM7782927
#   GSM7782928
#   GSM7782929
#   GSM7782930
#
# ================================================================


sample_groups <- data.frame(
  
  Sample_ID = c(
    
    "GSM7782916",
    "GSM7782917",
    
    "GSM7782918",
    "GSM7782919",
    "GSM7782920",
    "GSM7782921",
    "GSM7782922",
    
    "GSM7782923",
    "GSM7782924",
    "GSM7782925",
    "GSM7782926",
    "GSM7782927",
    "GSM7782928",
    "GSM7782929",
    "GSM7782930"
    
  ),
  
  Disease_Status = c(
    
    "Normal",
    "Normal",
    
    "Pathological_Aging",
    "Pathological_Aging",
    "Pathological_Aging",
    "Pathological_Aging",
    "Pathological_Aging",
    
    "Alzheimers_Disease",
    "Alzheimers_Disease",
    "Alzheimers_Disease",
    "Alzheimers_Disease",
    "Alzheimers_Disease",
    "Alzheimers_Disease",
    "Alzheimers_Disease",
    "Alzheimers_Disease"
    
  ),
  
  stringsAsFactors = FALSE
  
)


# ================================================================
# SECTION 10 — FIND SAMPLE ID COLUMN
# ================================================================
#
# The exact metadata column name may differ between H5AD files.
# Therefore we search for the column containing GSM IDs instead of
# assuming a particular column name.
#
# ================================================================

metadata_columns <- colnames(
  seurat_obj[[]]
)


gsm_pattern <- "^GSM77829(1[6-9]|2[0-9]|30)$"


gsm_hits <- sapply(
  
  metadata_columns,
  
  function(column_name) {
    
    values <- as.character(
      
      seurat_obj[[]][[column_name]]
      
    )
    
    any(
      grepl(
        gsm_pattern,
        values
      )
    )
    
  }
  
)


gsm_columns <- names(
  gsm_hits
)[gsm_hits]


if (length(gsm_columns) == 0) {
  
  stop(
    paste(
      "Could not find a metadata column containing",
      "GSM7782916-GSM7782930.",
      "\nRun:",
      "colnames(seurat_obj[[]])",
      "\nand inspect the metadata."
    )
  )
}


sample_column <- gsm_columns[1]


cat(
  "Detected sample column:",
  sample_column,
  "\n"
)


# ================================================================
# SECTION 11 — ASSIGN SAMPLE ID
# ================================================================

seurat_obj$Sample_ID <- as.character(
  
  seurat_obj[[]][[sample_column]]
  
)


# ================================================================
# SECTION 12 — ASSIGN DISEASE STATUS
# ================================================================

seurat_obj$Disease_Status <- sample_groups$Disease_Status[
  
  match(
    
    seurat_obj$Sample_ID,
    
    sample_groups$Sample_ID
    
  )
  
]


# Check disease groups

table(
  
  seurat_obj$Disease_Status,
  
  useNA = "ifany"
  
)


# Check samples

table(
  seurat_obj$Sample_ID
)


# Save initial object

saveRDS(
  
  seurat_obj,
  
  "data/processed/GSE243292_01_loaded.rds"
  
)


# ================================================================
# SECTION 13 — QUALITY CONTROL
# ================================================================
#
# Calculate mitochondrial percentage.
#
# Human mitochondrial genes generally begin with:
#
#   MT-
#
# ================================================================

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  
  seurat_obj,
  
  pattern = "^MT-"
  
)


# ================================================================
# SECTION 14 — QC VISUALIZATION
# ================================================================

qc_before <- VlnPlot(
  
  seurat_obj,
  
  features = c(
    
    "nFeature_RNA",
    
    "nCount_RNA",
    
    "percent.mt"
    
  ),
  
  ncol = 3,
  
  pt.size = 0.05
  
)


qc_before


ggsave(
  
  "figures/QC_before_filtering.png",
  
  qc_before,
  
  width = 12,
  
  height = 5
  
)


# ================================================================
# SECTION 15 — QC SCATTER PLOTS
# ================================================================

qc_scatter1 <- FeatureScatter(
  
  seurat_obj,
  
  feature1 = "nCount_RNA",
  
  feature2 = "nFeature_RNA"
  
)


qc_scatter2 <- FeatureScatter(
  
  seurat_obj,
  
  feature1 = "nCount_RNA",
  
  feature2 = "percent.mt"
  
)


qc_scatter <- qc_scatter1 + qc_scatter2


qc_scatter


ggsave(
  
  "figures/QC_scatter.png",
  
  qc_scatter,
  
  width = 12,
  
  height = 5
  
)


# ================================================================
# SECTION 16 — QC FILTERING
# ================================================================
#
# Baseline thresholds:
#
# nFeature_RNA >= 200
# nFeature_RNA <= 4000
# percent.mt < 5
#
# These should be checked against the QC plots and documented in
# the final report.
#
# ================================================================

cells_before_qc <- ncol(
  seurat_obj
)


seurat_obj <- subset(
  
  seurat_obj,
  
  subset =
    
    nFeature_RNA >= 200 &
    
    nFeature_RNA <= 4000 &
    
    percent.mt < 5
  
)


cells_after_qc <- ncol(
  seurat_obj
)


cat(
  "Cells before QC:",
  cells_before_qc,
  "\n"
)


cat(
  "Cells after QC:",
  cells_after_qc,
  "\n"
)


cat(
  "Retention:",
  round(
    100 * cells_after_qc /
      cells_before_qc,
    2
  ),
  "%\n"
)


# ================================================================
# SECTION 17 — QC AFTER FILTERING
# ================================================================

qc_after <- VlnPlot(
  
  seurat_obj,
  
  features = c(
    
    "nFeature_RNA",
    
    "nCount_RNA",
    
    "percent.mt"
    
  ),
  
  ncol = 3,
  
  pt.size = 0.05
  
)


qc_after


ggsave(
  
  "figures/QC_after_filtering.png",
  
  qc_after,
  
  width = 12,
  
  height = 5
  
)


# Save QC checkpoint

saveRDS(
  
  seurat_obj,
  
  "data/processed/GSE243292_02_QC.rds"
  
)
summary(seurat_obj$nFeature_RNA)
summary(seurat_obj$nCount_RNA)
summary(seurat_obj$percent.mt)
quantile(
  seurat_obj$nFeature_RNA,
  probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)
)


# ================================================================
# SECTION 18 — NORMALIZATION
# ================================================================
#
# LogNormalize is used for the baseline Seurat workflow.
#
# scale.factor = 10000
#
# ================================================================

seurat_obj <- NormalizeData(
  
  seurat_obj,
  
  normalization.method = "LogNormalize",
  
  scale.factor = 10000,
  
  verbose = FALSE
  
)


# ================================================================
# SECTION 19 — HIGHLY VARIABLE GENES
# ================================================================
#
# Use 2,000 highly variable genes.
#
# This keeps downstream memory use manageable.
#
# ================================================================

seurat_obj <- FindVariableFeatures(
  
  seurat_obj,
  
  selection.method = "vst",
  
  nfeatures = 2000,
  
  verbose = FALSE
  
)


# Plot variable genes

hvg_plot <- VariableFeaturePlot(
  seurat_obj
)


hvg_plot


ggsave(
  
  "figures/highly_variable_genes.png",
  
  hvg_plot,
  
  width = 8,
  
  height = 6
  
)


# Save HVG list

write.csv(
  
  data.frame(
    
    gene = VariableFeatures(
      seurat_obj
    )
    
  ),
  
  "results/highly_variable_genes.csv",
  
  row.names = FALSE
  
)


# ================================================================
# SECTION 20 — SCALING
# ================================================================
#
# IMPORTANT FOR 12 GB RAM:
#
# Scale only the 2,000 variable genes instead of all genes.
#
# ================================================================

variable_genes <- VariableFeatures(
  seurat_obj
)


seurat_obj <- ScaleData(
  
  seurat_obj,
  
  features = variable_genes,
  
  verbose = FALSE
  
)


# ================================================================
# SECTION 21 — PCA
# ================================================================
#
# Calculate 30 PCs.
#
# We initially use PCs 1-20 downstream.
#
# BUG FIX: RunPCA (irlba), FindClusters (Louvain), RunUMAP, and
# RunTSNE all involve stochastic steps. Without a fixed seed the
# cluster IDs, UMAP layout, and t-SNE layout are not reproducible
# between runs -- which matters here because cluster_annotation in
# SECTION 31 hardcodes cell-type labels to specific cluster numbers.
# set.seed() is called before each stochastic step below.
# ================================================================

set.seed(42)

seurat_obj <- RunPCA(
  
  seurat_obj,
  
  features = variable_genes,
  
  npcs = 30,
  
  verbose = FALSE
  
)


# PCA summary

print(
  seurat_obj[["pca"]],
  dims = 1:5,
  nfeatures = 10
)


# ================================================================
# SECTION 22 — ELBOW PLOT
# ================================================================

elbow <- ElbowPlot(
  
  seurat_obj,
  
  ndims = 30
  
)


elbow


ggsave(
  
  "figures/PCA_elbow_plot.png",
  
  elbow,
  
  width = 7,
  
  height = 5
  
)


# ================================================================
# SECTION 23 — CLUSTERING
# ================================================================
#
# Initial settings:
#
# PCs:
#   1:20
#
# Resolution:
#   0.5
#
# These can be adjusted after inspecting the biological structure.
#
# ================================================================
graphics.off()

dims_use <- 1:20

set.seed(42)

seurat_obj <- FindNeighbors(
  
  seurat_obj,
  
  reduction = "pca",
  
  dims = dims_use,
  
  verbose = FALSE
  
)

#names(seurat_obj@graphs)

seurat_obj <- FindClusters(
  
  seurat_obj,
  
  resolution = 0.5,
  
  verbose = FALSE
  
)


# Check cluster sizes

table(
  Idents(seurat_obj)
)


# ================================================================
# SECTION 24 — UMAP
# ================================================================

set.seed(42)

seurat_obj <- RunUMAP(
  
  seurat_obj,
  
  reduction = "pca",
  
  dims = dims_use,
  
  verbose = FALSE
  
)



umap_clusters <- DimPlot(
  
  seurat_obj,
  
  reduction = "umap",
  
  group.by = "seurat_clusters",
  
  label = TRUE,
  
  repel = TRUE
  
)


umap_clusters


ggsave(
  
  "figures/UMAP_clusters.png",
  
  umap_clusters,
  
  width = 10,
  
  height = 7
  
)


# ================================================================
# SECTION 25 — UMAP BY SAMPLE
# ================================================================
umap_sample <- DimPlot(
  
  seurat_obj,
  
  reduction = "umap",
  
  group.by = "Sample_ID"
  
)


umap_sample


ggsave(
  
  "figures/UMAP_by_sample.png",
  
  umap_sample,
  
  width = 10,
  
  height = 7
  
)


# ================================================================
# SECTION 26 — UMAP BY DISEASE STATUS
# ================================================================

umap_disease <- DimPlot(
  
  seurat_obj,
  
  reduction = "umap",
  
  group.by = "Disease_Status"
  
)


umap_disease


ggsave(
  
  "figures/UMAP_by_disease.png",
  
  umap_disease,
  
  width = 10,
  
  height = 7
  
)


table(
  seurat_obj$cell_type,
  seurat_obj$Disease_Status
)


# ================================================================
# SECTION 27 — t-SNE
# ================================================================

set.seed(42)

seurat_obj <- RunTSNE(
  
  seurat_obj,
  
  reduction = "pca",
  
  dims = dims_use,
  
  verbose = FALSE
  
)


tsne_plot <- DimPlot(
  
  seurat_obj,
  
  reduction = "tsne",
  
  group.by = "seurat_clusters",
  
  label = TRUE,
  
  repel = TRUE
  
)


tsne_plot


ggsave(
  
  "figures/tSNE_clusters.png",
  
  tsne_plot,
  
  width = 10,
  
  height = 7
  
)


# Save clustering checkpoint

saveRDS(
  
  seurat_obj,
  
  "data/processed/GSE243292_05_clustering.rds"
  
)


# ================================================================
# SECTION 28 — CELL-TYPE ANNOTATION
# ================================================================
#
# Canonical markers are used to identify likely cell types.
#
# IMPORTANT:
# The marker list is NOT an automatic annotation.
#
# You must inspect the DotPlot and FeaturePlots and determine
# which cluster corresponds to which cell type.
#
# ================================================================


cell_markers <- c(
  
  # Microglia
  "AIF1",
  "CX3CR1",
  "P2RY12",
  "C1QA",
  "C1QB",
  
  # Astrocytes
  "GFAP",
  "AQP4",
  "ALDH1L1",
  "SLC1A3",
  
  # Excitatory neurons
  "SLC17A7",
  "CAMK2A",
  "SATB2",
  
  # Inhibitory neurons
  "GAD1",
  "GAD2",
  "SLC6A1",
  
  # Oligodendrocytes
  "MBP",
  "PLP1",
  "MOG",
  "MOBP",
  
  # OPC
  "PDGFRA",
  "CSPG4",
  
  # Endothelial
  "CLDN5",
  "EMCN",
  "KDR"
  
)


# Keep genes actually present in dataset

cell_markers <- cell_markers[
  cell_markers %in%
    rownames(seurat_obj)
]


# ================================================================
# SECTION 29 — CELL MARKER DOTPLOT
# ================================================================

cell_marker_dotplot <- DotPlot(
  
  seurat_obj,
  
  features = cell_markers
  
) +
  
  RotatedAxis()


cell_marker_dotplot


ggsave(
  
  "figures/celltype_marker_DotPlot.png",
  
  cell_marker_dotplot,
  
  width = 12,
  
  height = 8
  
)


# ================================================================
# SECTION 30 — CELL MARKER FEATUREPLOTS
# ================================================================

cell_marker_features <- FeaturePlot(
  
  seurat_obj,
  
  features = cell_markers,
  
  min.cutoff = "q10",
  
  max.cutoff = "q90"
  
)


cell_marker_features


ggsave(
  
  "figures/celltype_FeaturePlots.png",
  
  cell_marker_features,
  
  width = 14,
  
  height = 10
  
)


# ============================================================
# ADDTIONAL STEP TO FIND MARKERS FOR ALL 26 CLUSTERS
# ============================================================

message("Finding cluster-specific markers...")

all_markers <- FindAllMarkers(
  
  object = seurat_obj,
  
  assay = "RNA",
  
  only.pos = TRUE,
  
  min.pct = 0.25,
  
  logfc.threshold = 0.25,
  
  return.thresh = 0.05,
  
  test.use = "wilcox",
  
  verbose = TRUE
  
)

message("Marker identification completed.")

exists("all_markers")

dim(all_markers)

head(all_markers)


# ============================================================
# TOP SIGNIFICANT MARKERS FOR EACH CLUSTER
# GSE243292
#
# Input:
#   all_markers
#
# Current dataset:
#   122,606 nuclei
#   26 clusters (0–25)
#
# IMPORTANT:
#   FindAllMarkers() has ALREADY been run.
#   This section only filters, ranks, extracts, and saves
#   the existing marker results.
# ============================================================


# ------------------------------------------------------------
# 1. Check the existing marker object
# ------------------------------------------------------------

dim(all_markers)

colnames(all_markers)

head(all_markers)


# ------------------------------------------------------------
# 2. Create results directory
# ------------------------------------------------------------

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)


# ------------------------------------------------------------
# 3. Identify the log-fold-change column
# ------------------------------------------------------------
#
# Seurat 5 normally uses:
#
#   avg_log2FC
#
# We check for it explicitly so that the code does not silently
# use the wrong column.
# ------------------------------------------------------------

if (!"avg_log2FC" %in% colnames(all_markers)) {
  
  stop(
    "Column 'avg_log2FC' was not found in all_markers."
  )
  
}


# ------------------------------------------------------------
# 4. Check the adjusted p-value column
# ------------------------------------------------------------

if (!"p_val_adj" %in% colnames(all_markers)) {
  
  stop(
    "Column 'p_val_adj' was not found in all_markers."
  )
  
}


# ------------------------------------------------------------
# 5. Check cluster column
# ------------------------------------------------------------

if (!"cluster" %in% colnames(all_markers)) {
  
  stop(
    "Column 'cluster' was not found in all_markers."
  )
  
}


# ------------------------------------------------------------
# 6. Filter significant positive markers
# ------------------------------------------------------------
#
# Criteria:
#
#   Adjusted P value <= 0.05
#   avg_log2FC >= 0.25
#
# These are reasonable starting criteria for annotation.
#
# Because all_markers already contains positive markers from
# FindAllMarkers(only.pos = TRUE), we do NOT need to rerun
# FindAllMarkers().
# ------------------------------------------------------------

significant_markers <- all_markers[
  
  all_markers$p_val_adj <= 0.05 &
    
    all_markers$avg_log2FC >= 0.25,
  
  ,
  drop = FALSE
  
]


# Check result

dim(significant_markers)


# ------------------------------------------------------------
# 7. Save complete significant marker table
# ------------------------------------------------------------

write.csv(
  
  significant_markers,
  
  "results/significant_cluster_markers.csv",
  
  row.names = FALSE
  
)


# ============================================================
# TOP 20 MARKERS PER CLUSTER
# ============================================================
#
# Ranking:
#
#   1. Adjusted P value
#   2. Average log2 fold-change
#
# The strongest markers should have:
#
#   low p_val_adj
#   high avg_log2FC
#
# ============================================================


top20_markers <- significant_markers |>
  
  dplyr::group_by(cluster) |>
  
  dplyr::arrange(
    p_val_adj,
    dplyr::desc(avg_log2FC),
    .by_group = TRUE
  ) |>
  
  dplyr::slice_head(n = 20) |>
  
  dplyr::ungroup()


# ------------------------------------------------------------
# Save top 20 markers
# ------------------------------------------------------------

write.csv(
  
  top20_markers,
  
  "results/top20_significant_markers_per_cluster.csv",
  
  row.names = FALSE
  
)


# ------------------------------------------------------------
# Check number of markers obtained per cluster
# ------------------------------------------------------------

table(
  top20_markers$cluster
)


# ============================================================
# TOP 10 MARKERS PER CLUSTER
# ============================================================
#
# This is a smaller annotation-friendly table.
# ============================================================


top10_markers <- top20_markers |>
  
  dplyr::group_by(cluster) |>
  
  dplyr::slice_head(n = 10) |>
  
  dplyr::ungroup()


# ------------------------------------------------------------
# Save top 10 markers
# ------------------------------------------------------------

write.csv(
  
  top10_markers,
  
  "results/top10_significant_markers_per_cluster.csv",
  
  row.names = FALSE
  
)


# ============================================================
# CREATE A SIMPLE ANNOTATION TABLE
# ============================================================
#
# This table is easier to inspect manually than the complete
# differential-expression output.
#
# It contains:
#
#   cluster
#   gene
#   avg_log2FC
#   pct.1
#   pct.2
#   p_val_adj
#
# ============================================================


annotation_marker_table <- top20_markers |>
  
  dplyr::select(
    cluster,
    gene,
    avg_log2FC,
    pct.1,
    pct.2,
    p_val_adj
  )


write.csv(
  
  annotation_marker_table,
  
  "results/annotation_marker_table_top20.csv",
  
  row.names = FALSE
  
)



# Create top 10 significant markers per cluster
top10_markers <- all_markers |>
  dplyr::filter(
    p_val_adj <= 0.05,
    avg_log2FC >= 0.25
  ) |>
  dplyr::group_by(cluster) |>
  dplyr::arrange(
    p_val_adj,
    dplyr::desc(avg_log2FC),
    .by_group = TRUE
  ) |>
  dplyr::slice_head(n = 10) |>
  dplyr::ungroup()


# Create wide annotation table
top10_wide <- top10_markers |>
  dplyr::group_by(cluster) |>
  dplyr::mutate(
    marker_rank = dplyr::row_number()
  ) |>
  dplyr::ungroup() |>
  dplyr::select(
    cluster,
    marker_rank,
    gene
  ) |>
  tidyr::pivot_wider(
    names_from = marker_rank,
    values_from = gene,
    names_prefix = "Marker_"
  )


# Save both tables
write.csv(
  top10_markers,
  "results/top10_significant_markers_per_cluster.csv",
  row.names = FALSE
)

write.csv(
  top10_wide,
  "results/top10_markers_annotation_wide.csv",
  row.names = FALSE
)





# ============================================================
# CREATE ONE ROW PER CLUSTER
# WITH TOP 10 GENES
# ============================================================
#
# This format is especially convenient for quickly reviewing
# the 26 clusters.
#
# Example:
#
# cluster | Marker1 | Marker2 | Marker3 | ...
#
# ============================================================


top10_wide <- top10_markers |>
  
  dplyr::group_by(cluster) |>
  
  dplyr::mutate(
    
    marker_rank = dplyr::row_number()
    
  ) |>
  
  dplyr::ungroup() |>
  
  dplyr::select(
    cluster,
    marker_rank,
    gene
  ) |>
  
  tidyr::pivot_wider(
    
    names_from = marker_rank,
    
    values_from = gene,
    
    names_prefix = "Marker_"
    
  )


# ------------------------------------------------------------
# Save wide annotation table
# ------------------------------------------------------------

write.csv(
  
  top10_wide,
  
  "results/top10_markers_annotation_wide.csv",
  
  row.names = FALSE
  
)


# ============================================================
# PRINT THE TOP 10 MARKERS FOR ALL 26 CLUSTERS
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "TOP MARKERS FOR EACH CLUSTER\n"
)

cat(
  "============================================\n"
)


for (cl in sort(
  unique(
    top10_markers$cluster
  )
)) {
  
  cat(
    "\nCluster ",
    cl,
    ":\n",
    sep = ""
  )
  
  genes <- top10_markers |>
    
    dplyr::filter(
      cluster == cl
    ) |>
    
    dplyr::pull(gene)
  
  cat(
    paste(
      genes,
      collapse = ", "
    ),
    "\n"
  )
  
}


# ============================================================
# FINAL SUMMARY
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "MARKER EXTRACTION COMPLETED\n"
)

cat(
  "============================================\n"
)

cat(
  "Original marker rows:",
  nrow(all_markers),
  "\n"
)

cat(
  "Significant marker rows:",
  nrow(significant_markers),
  "\n"
)

cat(
  "Clusters detected:",
  length(
    unique(
      all_markers$cluster
    )
  ),
  "\n"
)

cat(
  "\nFiles created:\n"
)

cat(
  "results/significant_cluster_markers.csv\n"
)

cat(
  "results/top20_significant_markers_per_cluster.csv\n"
)

cat(
  "results/top10_significant_markers_per_cluster.csv\n"
)

cat(
  "results/annotation_marker_table_top20.csv\n"
)

cat(
  "results/top10_markers_annotation_wide.csv\n"
)

cat(
  "\nUse these together with:\n"
)

cat(
  "figures/celltype_marker_DotPlot.png\n"
)

cat(
  "============================================\n"
)

# ================================================================
# SECTION 31 — CREATE CLUSTER ANNOTATION
# ================================================================
#
# IMPORTANT:
# Replace the example mapping below with your actual cluster
# identities after inspecting the marker plots.
#
# Example:
#
# cluster_annotation <- c(
#   "0" = "Microglia",
#   "1" = "Astrocytes",
#   "2" = "Oligodendrocytes"
# )
#
# ================================================================


# FIRST inspect:
#
# celltype_marker_DotPlot.png
#
# Then create your mapping.
#
# ---------------------------------------------------------------


cluster_annotation <- c(
  "0"  = "Oligodendrocytes",
  "1"  = "Oligodendrocytes",
  "2"  = "Astrocytes",
  "3"  = "Excitatory_Neurons",
  "4"  = "Excitatory_Neurons",
  "5"  = "OPC",
  "6"  = "Inhibitory_Neurons",
  "7"  = "Excitatory_Neurons",
  "8"  = "Excitatory_Neurons",
  "9"  = "Inhibitory_Neurons",
  "10" = "Excitatory_Neurons",
  "11" = "Microglia",
  "12" = "Excitatory_Neurons",
  "13" = "Unresolved",
  "14" = "Inhibitory_Neurons",
  "15" = "Astrocytes",
  "16" = "Inhibitory_Neurons",
  "17" = "Inhibitory_Neurons",
  "18" = "Endothelial",
  "19" = "Unresolved",
  "20" = "Astrocytes",
  "21" = "OPC",
  "22" = "Microglia",
  "23" = "OPC",
  "24" = "Inhibitory_Neurons",
  "25" = "Unresolved"
)


# Check whether annotation has been entered

if (length(cluster_annotation) == 0) {
  
  warning(
    paste(
      "No cluster annotation has been entered.",
      "\nInspect the marker plots and fill cluster_annotation."
    )
  )
  
} else {
  
  
  # Assign cell types
  
  seurat_obj$cell_type <- unname(
    
    cluster_annotation[
      as.character(
        seurat_obj$seurat_clusters
      )
    ]
    
  )
  
  
  # ==============================================================
  # ANNOTATED UMAP
  # ==============================================================
  
  annotated_umap <- DimPlot(
    
    seurat_obj,
    
    reduction = "umap",
    
    group.by = "cell_type",
    
    label = TRUE,
    
    repel = TRUE
    
  )
  
  
  annotated_umap
  
  
  ggsave(
    
    "figures/UMAP_annotated_cell_types.png",
    
    annotated_umap,
    
    width = 10,
    
    height = 7
    
  )
  
  
  # Cell-type counts
  
  celltype_counts <- as.data.frame(
    
    table(
      seurat_obj$cell_type
    )
    
  )
  
  
  colnames(
    celltype_counts
  ) <- c(
    
    "Cell_Type",
    
    "Nuclei"
    
  )
  
  
  write.csv(
    
    celltype_counts,
    
    "results/celltype_counts.csv",
    
    row.names = FALSE
    
  )
  
}


# ================================================================
# SECTION 32 — MARKER GENE IDENTIFICATION
# ================================================================
#
# FindAllMarkers identifies genes enriched in each cluster.
#
# ================================================================

markers <- FindAllMarkers(
  
  seurat_obj,
  
  only.pos = TRUE,
  
  min.pct = 0.25,
  
  logfc.threshold = 0.25,
  
  verbose = FALSE
  
)


# Save all markers

write.csv(
  
  markers,
  
  "results/all_cluster_markers.csv",
  
  row.names = FALSE
  
)


# ================================================================
# SECTION 33 — TOP 10 MARKERS PER CLUSTER
# ================================================================

top_markers <- markers %>%
  
  group_by(cluster) %>%
  
  slice_max(
    
    order_by = avg_log2FC,
    
    n = 10,
    
    with_ties = FALSE
    
  )


write.csv(
  
  top_markers,
  
  "results/top10_cluster_markers.csv",
  
  row.names = FALSE
  
)


# ================================================================
# SECTION 34 — MARKER HEATMAP
# ================================================================

heatmap_genes <- unique(

    top20_markers$gene

)

# Limit heatmap size for Avalible RAM

heatmap_genes <- head(

    heatmap_genes,
  
    50
)

marker_heatmap <- DoHeatmap(

    seurat_obj,
  
    features = heatmap_genes

) +
  NoLegend()


ggsave(

    "figures/top_marker_heatmap.png",
  
    marker_heatmap,
  
    width = 12,
  
    height = 10

)
#=================================================================
#saving checkpoint 
#=================================================================
saveRDS(
  seurat_obj,
  "data/processed/GSE243292_annotated_Seurat.rds"
)


# ================================================================
# DISEASE STATUS ANNOTATION — CONSISTENCY CHECK (BUG FIX)
# ================================================================
#
# BUG (fixed): this section used to *re-derive* Disease_Status from
# a second, independently hardcoded numeric map keyed on
# seurat_obj$sampleID (values "1"-"15"). That created two separate,
# unverified sources of truth for the same field:
#
#   1) Disease_Status assigned in SECTION 12 from sample_groups,
#      keyed on the GSM-based Sample_ID column.
#   2) Disease_Status re-assigned here from a second hardcoded
#      lookup keyed on sampleID.
#
# If sampleID's 1-15 ordering ever did not exactly match the GSM
# order assumed in SECTION 9, this would silently overwrite correct
# disease labels with incorrect ones, corrupting every downstream
# analysis (composition, DE, trajectory, CellChat).
#
# FIX: do not re-derive Disease_Status a second time. Instead,
# verify that the sampleID <-> Sample_ID <-> Disease_Status mapping
# already assigned in SECTION 12 is internally consistent (one
# Disease_Status per sampleID, no missing values). This preserves a
# single source of truth for Disease_Status.
# ================================================================

sample_lookup <- unique(
  seurat_obj[[]][, c("sampleID", "Sample_ID", "Disease_Status")]
)

sample_lookup <- sample_lookup[
  order(sample_lookup$sampleID),
]

print(sample_lookup)

if (anyNA(sample_lookup$Disease_Status)) {

  stop(
    paste(
      "Some sampleID/Sample_ID combinations have no assigned",
      "Disease_Status. Check the sample_groups table in SECTION 9."
    )
  )

}

if (any(table(sample_lookup$sampleID) > 1)) {

  stop(
    paste(
      "A single sampleID maps to more than one Sample_ID or",
      "Disease_Status. This indicates a metadata inconsistency",
      "between sampleID and the GSM-based Sample_ID."
    )
  )

}

cat(
  "Disease_Status verified consistent across",
  nrow(sample_lookup),
  "samples (single source of truth from SECTION 12).\n"
)


table(seurat_obj$Disease_Status)

table(
  seurat_obj$sampleID,
  seurat_obj$Disease_Status
)


#=================================================================
#Disease × Cell-type proportions 
#=================================================================

celltype_disease_prop <- prop.table(
  table(
    seurat_obj$cell_type,
    seurat_obj$Disease_Status
  ),
  margin = 2
)

round(
  100 * celltype_disease_prop,
  2
)


write.csv(
  round(100 * celltype_disease_prop, 2),
  "results/celltype_composition_by_disease.csv"
)

saveRDS(
  seurat_obj,
  "data/processed/GSE243292_disease_annotated.rds"
)


# ================================================================
# SECTION 35 — CELL-TYPE-SPECIFIC MARKERS
# ================================================================
#
# This section only runs if cell_type annotation exists.
#
# ================================================================

if ("cell_type" %in%
    colnames(seurat_obj[[]])) {
  
  
  Idents(seurat_obj) <- "cell_type"
  
  
  celltype_markers <- FindAllMarkers(
    
    seurat_obj,
    
    only.pos = TRUE,
    
    min.pct = 0.25,
    
    logfc.threshold = 0.25,
    
    verbose = FALSE
    
  )
  
  
  write.csv(
    
    celltype_markers,
    
    "results/all_celltype_markers.csv",
    
    row.names = FALSE
    
  )
  
  
  top_celltype_markers <- celltype_markers %>%
    
    group_by(cluster) %>%
    
    slice_max(
      
      order_by = avg_log2FC,
      
      n = 10,
      
      with_ties = FALSE
      
    )
  
  
  write.csv(
    
    top_celltype_markers,
    
    "results/top10_celltype_markers.csv",
    
    row.names = FALSE
    
  )
  
}

# ================================================================
# SECTION 36 — SAMPLE-AWARE PSEUDOBULK PREPARATION
# ================================================================
#
# Purpose:
# Aggregate raw RNA counts at the SAMPLE level within each
# annotated cell type.
#
# Biological replicates:
#   Normal              = samples 1–2
#   Pathological Aging  = samples 3–7
#   Alzheimer's Disease = samples 8–15
#
# This section prepares pseudobulk data.
# Statistical differential expression will be performed separately.
#
# ================================================================


# ================================================================
# SECTION 36.1 — CHECK SAMPLE / DISEASE MAPPING
# ================================================================

sample_info <- unique(
  seurat_obj[[]][,
                 c(
                   "sampleID",
                   "Disease_Status"
                 )
  ]
)

sample_info <- sample_info[
  order(sample_info$sampleID),
]

print(sample_info)


# ================================================================
# SECTION 36.2 — CHECK CELL-TYPE / SAMPLE COUNTS
# ================================================================

celltype_sample_table <- table(
  seurat_obj$cell_type,
  seurat_obj$sampleID
)

print(
  celltype_sample_table
)

write.csv(
  as.data.frame(
    celltype_sample_table
  ),
  "results/celltype_by_sample_counts.csv",
  row.names = FALSE
)


# ================================================================
# SECTION 36.3 — DEFINE CELL TYPES
# ================================================================
#
# Unresolved cells are excluded from the disease DE analysis.
#
# ================================================================

cell_types <- setdiff(
  unique(
    na.omit(
      seurat_obj$cell_type
    )
  ),
  "Unresolved"
)

print(
  cell_types
)


# ================================================================
# SECTION 36.4 — PREPARE RAW COUNTS
# ================================================================

metadata <- seurat_obj[[]]

counts_matrix <- GetAssayData(
  seurat_obj,
  assay = "RNA",
  layer = "counts"
)


# Storage object for pseudobulk results

pseudobulk_results <- list()


# ================================================================
# SECTION 36.5 — CREATE PSEUDOBULK DATA
# ================================================================

for (ct in cell_types) {
  
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "Processing cell type:",
    ct,
    "\n"
  )
  
  cat(
    "============================================\n"
  )
  
  
  # ------------------------------------------------
  # Get cells belonging to this cell type
  # ------------------------------------------------
  
  cells_ct <- rownames(
    metadata[
      !is.na(metadata$cell_type) &
        metadata$cell_type == ct,
      ,
      drop = FALSE
    ]
  )
  
  
  # Make sure cells exist in count matrix
  
  cells_ct <- intersect(
    cells_ct,
    colnames(counts_matrix)
  )
  
  
  if (
    length(cells_ct) == 0
  ) {
    
    cat(
      "No cells found. Skipping.\n"
    )
    
    next
    
  }
  
  
  # ------------------------------------------------
  # Metadata for current cell type
  # ------------------------------------------------
  
  meta_ct <- metadata[
    cells_ct,
    ,
    drop = FALSE
  ]
  
  
  # ------------------------------------------------
  # Samples represented in this cell type
  # ------------------------------------------------
  
  samples_ct <- sort(
    unique(
      meta_ct$sampleID
    )
  )
  
  
  # ------------------------------------------------
  # Aggregate counts by sample
  # ------------------------------------------------
  
  pb_list <- list()
  
  
  for (sid in samples_ct) {
    
    
    sample_cells <- rownames(
      meta_ct[
        meta_ct$sampleID == sid,
        ,
        drop = FALSE
      ]
    )
    
    
    sample_cells <- intersect(
      sample_cells,
      cells_ct
    )
    
    
    if (
      length(sample_cells) == 0
    ) {
      
      next
      
    }
    
    
    sample_counts <- counts_matrix[
      ,
      sample_cells,
      drop = FALSE
    ]
    
    
    pb_list[[as.character(sid)]] <-
      Matrix::rowSums(
        sample_counts
      )
    
  }
  
  
  # ------------------------------------------------
  # Combine sample-level counts
  # ------------------------------------------------
  
  pb_counts <- do.call(
    cbind,
    pb_list
  )
  
  
  # ------------------------------------------------
  # Create pseudobulk metadata
  # ------------------------------------------------
  
  pb_samples <- colnames(
    pb_counts
  )
  
  
  pb_metadata <- data.frame(
    
    sampleID = as.numeric(
      pb_samples
    ),
    
    row.names = pb_samples
    
  )
  
  
  # Add disease status
  
  pb_metadata$Disease_Status <-
    sample_info$Disease_Status[
      match(
        pb_metadata$sampleID,
        sample_info$sampleID
      )
    ]
  
  
  # ------------------------------------------------
  # Keep Alzheimer's Disease and Normal only
  # ------------------------------------------------
  
  keep_samples <-
    pb_metadata$Disease_Status %in%
    c(
      "Alzheimers_Disease",
      "Normal"
    )
  
  
  pb_counts <- pb_counts[
    ,
    keep_samples,
    drop = FALSE
  ]
  
  
  pb_metadata <- pb_metadata[
    keep_samples,
    ,
    drop = FALSE
  ]
  
  
  # ------------------------------------------------
  # Check biological replicates
  # ------------------------------------------------
  
  group_counts <- table(
    pb_metadata$Disease_Status
  )
  
  
  cat(
    "Sample counts:\n"
  )
  
  print(
    group_counts
  )
  
  
  # Need both groups
  
  if (
    !all(
      c(
        "Alzheimers_Disease",
        "Normal"
      ) %in%
      names(group_counts)
    )
  ) {
    
    cat(
      "Both disease groups are not available. Skipping.\n"
    )
    
    next
    
  }
  
  
  # Need at least two samples per group
  
  if (
    any(
      group_counts[
        c(
          "Alzheimers_Disease",
          "Normal"
        )
      ] < 2
    )
  ) {
    
    cat(
      "Insufficient biological replicates. Skipping.\n"
    )
    
    next
    
  }
  
  
  # ------------------------------------------------
  # Create safe filename
  # ------------------------------------------------
  
  safe_name <- gsub(
    "[^A-Za-z0-9_]+",
    "_",
    ct
  )
  
  
  # ------------------------------------------------
  # Save pseudobulk counts
  # ------------------------------------------------
  
  write.csv(
    as.matrix(
      pb_counts
    ),
    paste0(
      "results/pseudobulk_counts_",
      safe_name,
      ".csv"
    )
  )
  
  
  # ------------------------------------------------
  # Save pseudobulk metadata
  # ------------------------------------------------
  
  write.csv(
    pb_metadata,
    paste0(
      "results/pseudobulk_metadata_",
      safe_name,
      ".csv"
    ),
    row.names = TRUE
  )
  
  
  # ------------------------------------------------
  # Store in R
  # ------------------------------------------------
  
  pseudobulk_results[[ct]] <- list(
    
    counts = pb_counts,
    
    metadata = pb_metadata
    
  )
  
  
  cat(
    "Pseudobulk completed for:",
    ct,
    "\n"
  )
  
}


# ================================================================
# SECTION 36.6 — PSEUDOBULK SUMMARY
# ================================================================

cat(
  "\n============================================\n"
)

cat(
  "PSEUDOBULK PREPARATION COMPLETED\n"
)

cat(
  "============================================\n"
)

cat(
  "Cell types processed:",
  length(
    pseudobulk_results
  ),
  "\n"
)

cat(
  "\nResults saved under:\n"
)

cat(
  "results/\n"
)

cat(
  "\n============================================\n"
)


# ================================================================
# SECTION 37 — PSEUDOBULK DIFFERENTIAL EXPRESSION
# ================================================================
#
# Alzheimer's Disease vs Normal
#
# Uses sample-level pseudobulk counts.
#
# Biological replicates:
#   Normal              = 2 samples
#   Alzheimer's Disease = 8 samples
#
# Pathological Aging is not included in this comparison.
#
# ================================================================


# ------------------------------------------------
# SECTION 37.1 — LOAD edgeR
# ------------------------------------------------

if (!requireNamespace("edgeR", quietly = TRUE)) {
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  BiocManager::install("edgeR")
}

library(edgeR)


# ------------------------------------------------
# SECTION 37.2 — CREATE OUTPUT DIRECTORY
# ------------------------------------------------

dir.create(
  "results/pseudobulk_DE",
  showWarnings = FALSE,
  recursive = TRUE
)


# ------------------------------------------------
# SECTION 37.3 — INITIALIZE RESULTS
# ------------------------------------------------

pseudobulk_DE_results <- list()


# ------------------------------------------------
# SECTION 37.4 — RUN DE FOR EACH CELL TYPE
# ------------------------------------------------

for (ct in names(pseudobulk_results)) {
  
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "Running DE for:",
    ct,
    "\n"
  )
  
  cat(
    "============================================\n"
  )
  
  
  pb <- pseudobulk_results[[ct]]
  
  counts <- pb$counts
  
  metadata_pb <- pb$metadata
  
  
  # ------------------------------------------------
  # Disease status
  # ------------------------------------------------
  
  metadata_pb$Disease_Status <- factor(
    
    metadata_pb$Disease_Status,
    
    levels = c(
      "Normal",
      "Alzheimers_Disease"
    )
    
  )
  
  
  # ------------------------------------------------
  # Create edgeR object
  # ------------------------------------------------
  
  dge <- DGEList(
    
    counts = counts,
    
    group = metadata_pb$Disease_Status
    
  )
  
  
  # ------------------------------------------------
  # Filter low-expression genes
  # ------------------------------------------------
  
  keep_genes <- filterByExpr(
    
    dge,
    
    group = metadata_pb$Disease_Status
    
  )
  
  
  dge <- dge[
    
    keep_genes,
    
    ,
    
    keep.lib.sizes = FALSE
    
  ]
  
  
  cat(
    "Genes retained:",
    nrow(dge),
    "\n"
  )
  
  
  # ------------------------------------------------
  # Normalize
  # ------------------------------------------------
  
  dge <- calcNormFactors(
    dge
  )
  
  
  # ------------------------------------------------
  # Design matrix
  # ------------------------------------------------
  
  design <- model.matrix(
    
    ~ Disease_Status,
    
    data = metadata_pb
    
  )
  
  
  # ------------------------------------------------
  # Estimate dispersion
  # ------------------------------------------------
  
  dge <- estimateDisp(
    
    dge,
    
    design
    
  )
  
  
  # ------------------------------------------------
  # Fit model
  # ------------------------------------------------
  
  fit <- glmQLFit(
    
    dge,
    
    design,
    
    robust = TRUE
    
  )
  
  
  # ------------------------------------------------
  # Alzheimer's vs Normal
  # ------------------------------------------------
  
  qlf <- glmQLFTest(
    
    fit,
    
    coef = "Disease_StatusAlzheimers_Disease"
    
  )
  
  
  # ------------------------------------------------
  # Extract results
  # ------------------------------------------------
  
  deg <- topTags(
    
    qlf,
    
    n = Inf
    
  )$table
  
  
  deg$gene <- rownames(
    deg
  )
  
  
  deg <- deg[
    
    ,
    
    c(
      "gene",
      setdiff(
        colnames(deg),
        "gene"
      )
    )
    
  ]
  
  
  # ------------------------------------------------
  # Significance classification
  # ------------------------------------------------
  
  deg$significant <- ifelse(
    
    deg$FDR <= 0.05 &
      abs(deg$logFC) >= 0.25,
    
    "Significant",
    
    "Not_Significant"
    
  )
  
  
  # ------------------------------------------------
  # Safe filename
  # ------------------------------------------------
  
  safe_name <- gsub(
    
    "[^A-Za-z0-9_]+",
    
    "_",
    
    ct
    
  )
  
  
  # ------------------------------------------------
  # Save complete DE table
  # ------------------------------------------------
  
  write.csv(
    
    deg,
    
    paste0(
      "results/pseudobulk_DE/DE_",
      safe_name,
      "_Alzheimers_vs_Normal.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  # ------------------------------------------------
  # Save significant genes
  # ------------------------------------------------
  
  significant_deg <- deg[
    
    deg$FDR <= 0.05 &
      abs(deg$logFC) >= 0.25,
    
    ,
    
    drop = FALSE
    
  ]
  
  
  write.csv(
    
    significant_deg,
    
    paste0(
      "results/pseudobulk_DE/significant_DE_",
      safe_name,
      "_Alzheimers_vs_Normal.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  # ------------------------------------------------
  # Store results
  # ------------------------------------------------
  
  pseudobulk_DE_results[[ct]] <- deg
  
  
  cat(
    "Significant genes:",
    nrow(significant_deg),
    "\n"
  )
  
}


# ================================================================
# SECTION 37.5 — DE SUMMARY
# ================================================================

de_summary <- data.frame()


for (ct in names(pseudobulk_DE_results)) {
  
  
  deg <- pseudobulk_DE_results[[ct]]
  
  
  sig <- deg[
    
    deg$FDR <= 0.05 &
      abs(deg$logFC) >= 0.25,
    
    ,
    
    drop = FALSE
    
  ]
  
  
  de_summary <- rbind(
    
    de_summary,
    
    data.frame(
      
      Cell_Type = ct,
      
      Total_Genes = nrow(deg),
      
      Significant = nrow(sig),
      
      Up_in_AD = sum(
        sig$logFC > 0
      ),
      
      Down_in_AD = sum(
        sig$logFC < 0
      )
      
    )
    
  )
  
}


write.csv(
  
  de_summary,
  
  "results/pseudobulk_DE/DE_summary.csv",
  
  row.names = FALSE
  
)


print(
  de_summary
)


# ================================================================
# SECTION 37.6 — COMPLETION
# ================================================================

cat(
  "\n============================================\n"
)

cat(
  "PSEUDOBULK DIFFERENTIAL EXPRESSION COMPLETED\n"
)

cat(
  "Comparison: Alzheimer's Disease vs Normal\n"
)

cat(
  "Results: results/pseudobulk_DE/\n"
)

cat(
  "============================================\n"
)


# ================================================================
# SECTION 38 — SAVE CURRENT ANALYSIS OBJECT
# ================================================================

saveRDS(
  seurat_obj,
  "data/processed/GSE243292_disease_DE_Seurat.rds"
)

write.csv(
  seurat_obj[[]],
  "results/disease_DE_cell_metadata.csv",
  row.names = TRUE
)



# ================================================================
# SECTION 39 — ANALYSIS SUMMARY
# ================================================================

cat(
  "\n\n============================================\n"
)

cat(
  "GSE243292 DISEASE ANALYSIS CHECKPOINT\n"
)

cat(
  "============================================\n"
)

cat(
  "Cells:",
  ncol(seurat_obj),
  "\n"
)

cat(
  "Genes:",
  nrow(seurat_obj),
  "\n"
)

cat(
  "Clusters:",
  length(
    unique(
      seurat_obj$seurat_clusters
    )
  ),
  "\n"
)

cat(
  "Cell types:",
  length(
    unique(
      na.omit(
        seurat_obj$cell_type
      )
    )
  ),
  "\n"
)

cat(
  "Disease groups:\n"
)

print(
  table(
    seurat_obj$Disease_Status
  )
)

cat(
  "\nPseudobulk DE results:\n",
  "results/pseudobulk_DE/\n"
)

cat(
  "\nSaved object:\n",
  "data/processed/GSE243292_disease_DE_Seurat.rds\n"
)

cat(
  "\n============================================\n"
)


# ================================================================
# SECTION 40 — TOP DIFFERENTIALLY EXPRESSED GENES
# ================================================================

for (ct in names(pseudobulk_DE_results)) {
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "CELL TYPE:",
    ct,
    "\n"
  )
  
  cat(
    "============================================\n"
  )
  
  
  deg <- pseudobulk_DE_results[[ct]]
  
  
  # Significant genes
  sig <- deg[
    deg$FDR <= 0.05 &
      abs(deg$logFC) >= 0.25,
    ,
    drop = FALSE
  ]
  
  
  # Sort by FDR
  sig <- sig[
    order(
      sig$FDR,
      -abs(sig$logFC)
    ),
    ,
    drop = FALSE
  ]
  
  
  print(
    head(
      sig,
      20
    )
  )
  
}



# ================================================================
# SECTION 41 — PSEUDOBULK DE VOLCANO PLOTS
# ================================================================

dir.create(
  "figures/pseudobulk_DE",
  showWarnings = FALSE,
  recursive = TRUE
)


for (ct in names(pseudobulk_DE_results)) {
  
  
  cat(
    "\nCreating volcano plot for:",
    ct,
    "\n"
  )
  
  
  deg <- pseudobulk_DE_results[[ct]]
  
  
  # ------------------------------------------------
  # Calculate -log10(FDR)
  # ------------------------------------------------
  
  deg$neg_log10_FDR <- -log10(
    pmax(
      deg$FDR,
      .Machine$double.xmin
    )
  )
  
  
  # ------------------------------------------------
  # Classify genes
  # ------------------------------------------------
  
  deg$Status <- "Not Significant"
  
  
  deg$Status[
    deg$FDR <= 0.05 &
      deg$logFC >= 0.25
  ] <- "Up in AD"
  
  
  deg$Status[
    deg$FDR <= 0.05 &
      deg$logFC <= -0.25
  ] <- "Down in AD"
  
  
  # ------------------------------------------------
  # Volcano plot
  # ------------------------------------------------
  
  volcano_plot <- ggplot2::ggplot(
    
    deg,
    
    ggplot2::aes(
      x = logFC,
      y = neg_log10_FDR
    )
    
  ) +
    
    ggplot2::geom_point(
      alpha = 0.6,
      size = 1
    ) +
    
    ggplot2::geom_vline(
      xintercept = c(
        -0.25,
        0.25
      ),
      linetype = "dashed"
    ) +
    
    ggplot2::geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    
    ggplot2::labs(
      
      title = paste(
        ct,
        "— Alzheimer's Disease vs Normal"
      ),
      
      x = "log2 Fold Change",
      
      y = "-log10(FDR)"
      
    ) +
    
    ggplot2::theme_classic()
  
  
  # ------------------------------------------------
  # Display
  # ------------------------------------------------
  
  print(
    volcano_plot
  )
  
  
  # ------------------------------------------------
  # Save
  # ------------------------------------------------
  
  safe_name <- gsub(
    "[^A-Za-z0-9_]+",
    "_",
    ct
  )
  
  
  ggplot2::ggsave(
    
    paste0(
      "figures/pseudobulk_DE/volcano_",
      safe_name,
      ".png"
    ),
    
    volcano_plot,
    
    width = 8,
    
    height = 6
  )
  
  
}


cat(
  "\n============================================\n"
)

cat(
  "SECTION 41 COMPLETED\n"
)

cat(
  "Volcano plots saved in:\n",
  "figures/pseudobulk_DE/\n"
)

cat(
  "============================================\n"
)


# ================================================================
# SECTION 42 — PSEUDOBULK DE HEATMAP
# ================================================================

# ------------------------------------------------
# Create output directory
# ------------------------------------------------

dir.create(
  "figures/pseudobulk_DE",
  showWarnings = FALSE,
  recursive = TRUE
)


# ------------------------------------------------
# Collect significant DE genes
# ------------------------------------------------

de_heatmap_genes <- c()


for (ct in names(pseudobulk_DE_results)) {
  
  
  deg <- pseudobulk_DE_results[[ct]]
  
  
  sig <- deg[
    deg$FDR <= 0.05 &
      abs(deg$logFC) >= 0.25,
    ,
    drop = FALSE
  ]
  
  
  # Take maximum 10 genes per cell type
  
  sig <- sig[
    order(
      sig$FDR,
      -abs(sig$logFC)
    ),
    ,
    drop = FALSE
  ]
  
  
  sig <- head(
    sig,
    10
  )
  
  
  de_heatmap_genes <- c(
    de_heatmap_genes,
    sig$gene
  )
  
}


# ------------------------------------------------
# Remove duplicate genes
# ------------------------------------------------

de_heatmap_genes <- unique(
  de_heatmap_genes
)


# ------------------------------------------------
# Keep only genes present in ScaleData
# ------------------------------------------------

scale_genes <- rownames(
  LayerData(
    seurat_obj,
    assay = "RNA",
    layer = "scale.data"
  )
)


de_heatmap_genes <- intersect(
  de_heatmap_genes,
  scale_genes
)


# ------------------------------------------------
# Limit number of genes for RAM
# ------------------------------------------------

de_heatmap_genes <- head(
  de_heatmap_genes,
  50
)


cat(
  "Genes used for heatmap:",
  length(de_heatmap_genes),
  "\n"
)


# ------------------------------------------------
# Check that genes exist
# ------------------------------------------------

if (length(de_heatmap_genes) == 0) {
  
  stop(
    "No DE genes were found in the RNA scale.data layer."
  )
}


# ------------------------------------------------
# Create heatmap
# ------------------------------------------------

de_heatmap <- DoHeatmap(
  
  seurat_obj,
  
  features = de_heatmap_genes,
  
  group.by = "Disease_Status"
  
) +
  
  NoLegend()


# ------------------------------------------------
# Display
# ------------------------------------------------

de_heatmap


# ------------------------------------------------
# Save
# ------------------------------------------------

ggsave(
  
  "figures/pseudobulk_DE/DE_heatmap_by_disease.png",
  
  de_heatmap,
  
  width = 12,
  
  height = 10
)


cat(
  "\n============================================\n"
)

cat(
  "SECTION 42 COMPLETED\n"
)

cat(
  "DE heatmap saved:\n",
  "figures/pseudobulk_DE/DE_heatmap_by_disease.png\n"
)

cat(
  "============================================\n"
)


# ================================================================
# SECTION 43 — GO BIOLOGICAL PROCESS ENRICHMENT
# ================================================================

# ------------------------------------------------
# Install/load required packages
# ------------------------------------------------

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

BiocManager::install(
  "enrichplot",
  ask = FALSE,
  update = FALSE
)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)


# ------------------------------------------------
# Create output directory
# ------------------------------------------------

dir.create(
  "results/pathway_enrichment",
  showWarnings = FALSE,
  recursive = TRUE
)


# ------------------------------------------------
# Initialize results
# ------------------------------------------------

GO_results <- list()


# ================================================================
# SECTION 43.1 — ENRICHMENT FOR EACH CELL TYPE
# ================================================================

for (ct in names(pseudobulk_DE_results)) {
  
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "PATHWAY ANALYSIS:",
    ct,
    "\n"
  )
  
  cat(
    "============================================\n"
  )
  
  
  deg <- pseudobulk_DE_results[[ct]]
  
  
  # ------------------------------------------------
  # Significant AD-upregulated genes
  # ------------------------------------------------
  
  up_genes <- deg$gene[
    
    deg$FDR <= 0.05 &
      deg$logFC > 0.25
    
  ]
  
  
  # ------------------------------------------------
  # Significant AD-downregulated genes
  # ------------------------------------------------
  
  down_genes <- deg$gene[
    
    deg$FDR <= 0.05 &
      deg$logFC < -0.25
    
  ]
  
  
  cat(
    "AD-upregulated genes:",
    length(up_genes),
    "\n"
  )
  
  cat(
    "AD-downregulated genes:",
    length(down_genes),
    "\n"
  )
  
  
  # ------------------------------------------------
  # Convert gene symbols to Entrez IDs
  # ------------------------------------------------
  
  up_entrez <- bitr(
    
    up_genes,
    
    fromType = "SYMBOL",
    
    toType = "ENTREZID",
    
    OrgDb = org.Hs.eg.db
    
  )
  
  
  down_entrez <- bitr(
    
    down_genes,
    
    fromType = "SYMBOL",
    
    toType = "ENTREZID",
    
    OrgDb = org.Hs.eg.db
    
  )
  
  
  # ------------------------------------------------
  # GO — AD UPREGULATED
  # ------------------------------------------------
  
  if (nrow(up_entrez) >= 3) {
    
    
    GO_up <- enrichGO(
      
      gene = unique(
        up_entrez$ENTREZID
      ),
      
      OrgDb = org.Hs.eg.db,
      
      keyType = "ENTREZID",
      
      ont = "BP",
      
      pAdjustMethod = "BH",
      
      pvalueCutoff = 0.05,
      
      qvalueCutoff = 0.2,
      
      readable = TRUE
      
    )
    
    
    if (!is.null(GO_up) &&
        nrow(as.data.frame(GO_up)) > 0) {
      
      write.csv(
        
        as.data.frame(GO_up),
        
        paste0(
          "results/pathway_enrichment/",
          gsub(
            "[^A-Za-z0-9_]+",
            "_",
            ct
          ),
          "_GO_BP_AD_up.csv"
        ),
        
        row.names = FALSE
        
      )
      
    }
    
  } else {
    
    cat(
      "Too few AD-upregulated genes for GO enrichment.\n"
    )
    
    GO_up <- NULL
    
  }
  
  
  # ------------------------------------------------
  # GO — AD DOWNREGULATED
  # ------------------------------------------------
  
  if (nrow(down_entrez) >= 3) {
    
    
    GO_down <- enrichGO(
      
      gene = unique(
        down_entrez$ENTREZID
      ),
      
      OrgDb = org.Hs.eg.db,
      
      keyType = "ENTREZID",
      
      ont = "BP",
      
      pAdjustMethod = "BH",
      
      pvalueCutoff = 0.05,
      
      qvalueCutoff = 0.2,
      
      readable = TRUE
      
    )
    
    
    if (!is.null(GO_down) &&
        nrow(as.data.frame(GO_down)) > 0) {
      
      write.csv(
        
        as.data.frame(GO_down),
        
        paste0(
          "results/pathway_enrichment/",
          gsub(
            "[^A-Za-z0-9_]+",
            "_",
            ct
          ),
          "_GO_BP_AD_down.csv"
        ),
        
        row.names = FALSE
        
      )
      
    }
    
  } else {
    
    cat(
      "Too few AD-downregulated genes for GO enrichment.\n"
    )
    
    GO_down <- NULL
    
  }
  
  
  # ------------------------------------------------
  # Store results
  # ------------------------------------------------
  
  GO_results[[ct]] <- list(
    
    up = GO_up,
    
    down = GO_down
    
  )
  
  
  cat(
    "GO enrichment completed.\n"
  )
  
}


# ================================================================
# SECTION 43.2 — COMPLETION
# ================================================================

cat(
  "\n============================================\n"
)

cat(
  "SECTION 43 COMPLETED\n"
)

cat(
  "GO Biological Process results saved in:\n"
)

cat(
  "results/pathway_enrichment/\n"
)

cat(
  "============================================\n"
)


# ================================================================
# SECTION 44 — GO ENRICHMENT VISUALIZATION
# ================================================================

# ------------------------------------------------
# Check available GO results
# ------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "GO ENRICHMENT RESULTS\n"
)

cat(
  "============================================\n"
)


for (ct in names(GO_results)) {
  
  cat(
    "\nCell type:",
    ct,
    "\n"
  )
  
  
  # AD-upregulated
  
  if (!is.null(GO_results[[ct]]$up)) {
    
    up_df <- as.data.frame(
      GO_results[[ct]]$up
    )
    
    cat(
      "AD-up GO terms:",
      nrow(up_df),
      "\n"
    )
    
  } else {
    
    cat(
      "AD-up GO terms: 0\n"
    )
    
  }
  
  
  # AD-downregulated
  
  if (!is.null(GO_results[[ct]]$down)) {
    
    down_df <- as.data.frame(
      GO_results[[ct]]$down
    )
    
    cat(
      "AD-down GO terms:",
      nrow(down_df),
      "\n"
    )
    
  } else {
    
    cat(
      "AD-down GO terms: 0\n"
    )
    
  }
  
}


# ================================================================
# SECTION 44.1 — PLOT AVAILABLE RESULTS
# ================================================================

dir.create(
  "figures/pathway_enrichment",
  showWarnings = FALSE,
  recursive = TRUE
)


for (ct in names(GO_results)) {
  
  
  safe_name <- gsub(
    "[^A-Za-z0-9_]+",
    "_",
    ct
  )
  
  
  # ------------------------------------------------
  # AD-upregulated GO
  # ------------------------------------------------
  
  if (!is.null(GO_results[[ct]]$up)) {
    
    up_df <- as.data.frame(
      GO_results[[ct]]$up
    )
    
    
    if (nrow(up_df) > 0) {
      
      p_up <- dotplot(
        GO_results[[ct]]$up,
        showCategory = min(
          15,
          nrow(up_df)
        )
      ) +
        ggplot2::ggtitle(
          paste(
            ct,
            "— GO BP: AD Upregulated"
          )
        )
      
      
      print(p_up)
      
      
      ggsave(
        paste0(
          "figures/pathway_enrichment/",
          safe_name,
          "_GO_BP_AD_up.png"
        ),
        p_up,
        width = 10,
        height = 7
      )
      
    }
    
  }
  
  
  # ------------------------------------------------
  # AD-downregulated GO
  # ------------------------------------------------
  
  if (!is.null(GO_results[[ct]]$down)) {
    
    down_df <- as.data.frame(
      GO_results[[ct]]$down
    )
    
    
    if (nrow(down_df) > 0) {
      
      p_down <- dotplot(
        GO_results[[ct]]$down,
        showCategory = min(
          15,
          nrow(down_df)
        )
      ) +
        ggplot2::ggtitle(
          paste(
            ct,
            "— GO BP: AD Downregulated"
          )
        )
      
      
      print(p_down)
      
      
      ggsave(
        paste0(
          "figures/pathway_enrichment/",
          safe_name,
          "_GO_BP_AD_down.png"
        ),
        p_down,
        width = 10,
        height = 7
      )
      
    }
    
  }
  
}


cat(
  "\n============================================\n"
)

cat(
  "SECTION 44 COMPLETED\n"
)

cat(
  "Pathway figures saved in:\n",
  "figures/pathway_enrichment/\n"
)

cat(
  "============================================\n"
)



# ================================================================
# SECTION 45 — INSTALL MONOCLE3
# ================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

BiocManager::install(
  c(
    "SingleCellExperiment",
    "batchelor",
    "HDF5Array",
    "ggrastr"
  ),
  ask = FALSE,
  update = FALSE
)

devtools::install_github(
  "cole-trapnell-lab/monocle3@v1.4.27",
  upgrade = "never"
)

library(monocle3)

cat("\nSECTION 45 COMPLETED — Monocle3 installed.\n")


# ================================================================
# SECTION 46 — PREPARE OPC → OLIGODENDROCYTE TRAJECTORY
# ================================================================

trajectory_cells <- colnames(seurat_obj)[
  seurat_obj$cell_type %in%
    c("OPC", "Oligodendrocytes")
]

cat(
  "Trajectory cells:",
  length(trajectory_cells),
  "\n"
)

print(
  table(
    seurat_obj$cell_type[trajectory_cells]
  )
)

trajectory_seurat <- subset(
  seurat_obj,
  cells = trajectory_cells
)

counts <- GetAssayData(
  trajectory_seurat,
  assay = "RNA",
  layer = "counts"
)

cell_metadata <- trajectory_seurat[[]]

gene_metadata <- data.frame(
  gene_short_name = rownames(counts),
  row.names = rownames(counts)
)

trajectory_cds <- new_cell_data_set(
  counts,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

cat("\nMonocle3 CDS created successfully.\n")
cat("Cells:", ncol(trajectory_cds), "\n")
cat("Genes:", nrow(trajectory_cds), "\n")


# ================================================================
# SECTION 47 — MONOCLE3 PREPROCESSING AND TRAJECTORY GRAPH
# ================================================================

trajectory_cds <- preprocess_cds(
  trajectory_cds,
  num_dim = 30
)

trajectory_cds <- reduce_dimension(
  trajectory_cds,
  reduction_method = "UMAP",
  preprocess_method = "PCA"
)

trajectory_cds <- cluster_cells(
  trajectory_cds,
  reduction_method = "UMAP"
)

cat("\nMonocle3 UMAP clustering completed.\n")

print(
  table(
    clusters(trajectory_cds)
  )
)

trajectory_umap <- plot_cells(
  trajectory_cds,
  color_cells_by = "cell_type",
  label_cell_groups = TRUE,
  label_leaves = TRUE,
  label_branch_points = TRUE
)

ggsave(
  "figures/trajectory_OPC_Oligodendrocytes_UMAP.png",
  trajectory_umap,
  width = 10,
  height = 7
)

trajectory_cds <- learn_graph(
  trajectory_cds,
  use_partition = TRUE
)

trajectory_graph <- plot_cells(
  trajectory_cds,
  color_cells_by = "cell_type",
  label_cell_groups = TRUE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  show_trajectory_graph = TRUE
)

ggsave(
  "figures/trajectory_OPC_Oligodendrocytes_graph.png",
  trajectory_graph,
  width = 10,
  height = 7
)

cat("\nSECTION 47 COMPLETED\n")


# ================================================================
# SECTION 48 — OPC → OLIGODENDROCYTE PSEUDOTIME
# ================================================================

opc_cells <- colnames(trajectory_cds)[
  colData(trajectory_cds)$cell_type == "OPC"
]

cat(
  "OPC root cells:",
  length(opc_cells),
  "\n"
)

trajectory_cds <- order_cells(
  trajectory_cds,
  root_cells = opc_cells
)

pseudotime_plot <- plot_cells(
  trajectory_cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  show_trajectory_graph = TRUE
)

ggsave(
  "figures/OPC_Oligodendrocyte_pseudotime.png",
  pseudotime_plot,
  width = 10,
  height = 7
)


# ================================================================
# SECTION 49 — VALIDATE PSEUDOTIME
# ================================================================

pseudotime_values <- pseudotime(
  trajectory_cds
)

pseudotime_summary <- data.frame(
  cell_type = colData(trajectory_cds)$cell_type,
  pseudotime = pseudotime_values
)

print(
  aggregate(
    pseudotime ~ cell_type,
    data = pseudotime_summary,
    FUN = median
  )
)

pseudotime_boxplot <- ggplot(
  pseudotime_summary,
  aes(
    x = cell_type,
    y = pseudotime
  )
) +
  geom_boxplot() +
  theme_classic() +
  labs(
    title = "OPC → Oligodendrocyte Pseudotime",
    x = "Cell type",
    y = "Pseudotime"
  )

ggsave(
  "figures/OPC_Oligodendrocyte_pseudotime_boxplot.png",
  pseudotime_boxplot,
  width = 8,
  height = 6
)

write.csv(
  pseudotime_summary,
  "results/OPC_Oligodendrocyte_pseudotime_values.csv",
  row.names = FALSE
)

cat("\nSECTION 49 COMPLETED\n")


# ================================================================
# SECTION 50 — TRAJECTORY-ASSOCIATED GENES
# ================================================================

gene_fits <- graph_test(
  trajectory_cds,
  neighbor_graph = "principal_graph",
  cores = 2
)

gene_fits <- gene_fits[
  order(gene_fits$q_value),
  ,
  drop = FALSE
]

significant_trajectory_genes <- gene_fits[
  gene_fits$q_value <= 0.05,
  ,
  drop = FALSE
]

write.csv(
  gene_fits,
  "results/trajectory_gene_test.csv",
  row.names = TRUE
)

write.csv(
  significant_trajectory_genes,
  "results/significant_trajectory_genes.csv",
  row.names = TRUE
)

cat(
  "\nSignificant trajectory genes:",
  nrow(significant_trajectory_genes),
  "\n"
)

cat("\nSECTION 50 COMPLETED\n")


# ================================================================
# SECTION 51 — TOP TRAJECTORY GENES AND VISUALIZATION
# ================================================================

top_trajectory_genes <- significant_trajectory_genes[
  order(significant_trajectory_genes$q_value),
  ,
  drop = FALSE
]

top_trajectory_genes <- head(
  top_trajectory_genes,
  50
)

write.csv(
  top_trajectory_genes,
  "results/top50_trajectory_genes.csv",
  row.names = TRUE
)

top_genes <- unique(
  top_trajectory_genes$gene_short_name[
    !is.na(top_trajectory_genes$gene_short_name) &
      top_trajectory_genes$gene_short_name != ""
  ]
)

top_genes <- head(
  top_genes,
  10
)

top_genes <- intersect(
  top_genes,
  rowData(trajectory_cds)$gene_short_name
)

trajectory_gene_cds <- trajectory_cds[
  rowData(trajectory_cds)$gene_short_name %in%
    top_genes
]

trajectory_gene_plot <- plot_genes_in_pseudotime(
  trajectory_gene_cds,
  min_expr = 0.5,
  ncol = 2
)

ggsave(
  "figures/top_trajectory_genes_pseudotime.png",
  trajectory_gene_plot,
  width = 12,
  height = 10
)

cat("\nSECTION 51 COMPLETED\n")


# ================================================================
# SECTION 52 — SAVE TRAJECTORY RESULTS
# ================================================================

dir.create(
  "results/trajectory",
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  pseudotime_summary,
  "results/trajectory/pseudotime_values.csv",
  row.names = FALSE
)

write.csv(
  top_trajectory_genes,
  "results/trajectory/top50_trajectory_genes.csv",
  row.names = TRUE
)

saveRDS(
  trajectory_cds,
  "data/processed/GSE243292_OPC_Oligodendrocyte_trajectory.rds"
)

cat(
  "\nSECTION 52 COMPLETED\n"
)

# ================================================================
# SECTION 53 — INSTALL CELLCHAT DEPENDENCIES
# ================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(
  "ComplexHeatmap",
  ask = FALSE,
  update = FALSE
)

library(ComplexHeatmap)

cat(
  "\nComplexHeatmap installed successfully.\n"
)


# ================================================================
# SECTION 54 — INSTALL CELLCHAT
# ================================================================

remotes::install_github(
  "jinworks/CellChat",
  upgrade = "never"
)

library(CellChat)

cat(
  "\nCellChat installed successfully.\n"
)

packageVersion("CellChat")

# ================================================================
# SECTION 55 — PREPARE CELLCHAT DATA
# ================================================================

# Remove unresolved cells

cellchat_cells <- rownames(
  seurat_obj[[]]
)[
  !is.na(seurat_obj$cell_type) &
    seurat_obj$cell_type != "Unresolved"
]

seurat_cellchat <- subset(
  seurat_obj,
  cells = cellchat_cells
)

# Separate disease groups

seurat_AD <- subset(
  seurat_cellchat,
  subset = Disease_Status == "Alzheimers_Disease"
)

seurat_Normal <- subset(
  seurat_cellchat,
  subset = Disease_Status == "Normal"
)

cat(
  "\nAD cells:",
  ncol(seurat_AD),
  "\n"
)

cat(
  "Normal cells:",
  ncol(seurat_Normal),
  "\n"
)

cat(
  "\nAD cell types:\n"
)

print(
  table(seurat_AD$cell_type)
)

cat(
  "\nNormal cell types:\n"
)

print(
  table(seurat_Normal$cell_type)
)

cat(
  "\nSECTION 55 COMPLETED\n"
)


# ================================================================
# SECTION 56 — CREATE CELLCHAT OBJECTS
# ================================================================

data.input.AD <- GetAssayData(
  seurat_AD,
  assay = "RNA",
  layer = "data"
)

data.input.Normal <- GetAssayData(
  seurat_Normal,
  assay = "RNA",
  layer = "data"
)

meta.AD <- seurat_AD[[]]
meta.Normal <- seurat_Normal[[]]

cellchat_AD <- createCellChat(
  object = data.input.AD,
  meta = meta.AD,
  group.by = "cell_type"
)

cellchat_Normal <- createCellChat(
  object = data.input.Normal,
  meta = meta.Normal,
  group.by = "cell_type"
)

cat(
  "\nCellChat objects created successfully.\n"
)

cat(
  "AD cells:",
  ncol(cellchat_AD@data),
  "\n"
)

cat(
  "Normal cells:",
  ncol(cellchat_Normal@data),
  "\n"
)

cat(
  "\nSECTION 56 COMPLETED\n"
)



# ================================================================
# SECTION 57 — LOAD HUMAN CELLCHAT DATABASE
# ================================================================

CellChatDB.human <- CellChatDB.human

cellchat_AD@DB <- CellChatDB.human
cellchat_Normal@DB <- CellChatDB.human

cat(
  "\nHuman CellChat database loaded.\n"
)
cellchat_AD <- subsetData(cellchat_AD)
cellchat_Normal <- subsetData(cellchat_Normal)

cat(
  "\nCellChat objects recreated successfully.\n"
)
cat(
  "Signaling interactions:",
  nrow(CellChatDB.human$interaction),
  "\n"
)

cat(
  "\nSECTION 57 COMPLETED\n"
)

# ================================================================
# SECTION 58 — INSTALL PRESTO
# ================================================================

remotes::install_github(
  "immunogenomics/presto",
  upgrade = "never"
)

library(presto)

cat(
  "\nPresto installed successfully.\n"
)


# ================================================================
# SECTION 59 — CELLCHAT OVEREXPRESSED GENES
# ================================================================

cellchat_AD <- identifyOverExpressedGenes(
  cellchat_AD
)

cellchat_AD <- identifyOverExpressedInteractions(
  cellchat_AD
)

cellchat_Normal <- identifyOverExpressedGenes(
  cellchat_Normal
)

cellchat_Normal <- identifyOverExpressedInteractions(
  cellchat_Normal
)

cat(
  "\nSECTION 59 COMPLETED\n"
)


# ================================================================
# SECTION 60 — COMPUTE COMMUNICATION PROBABILITIES
# ================================================================

cellchat_AD <- computeCommunProb(
  cellchat_AD
)

cellchat_Normal <- computeCommunProb(
  cellchat_Normal
)

cat(
  "\nCommunication probabi
  lities calculated.\n"
)

cat(
  "\nSECTION 60 COMPLETED\n"
)


# ================================================================
# SECTION 61 — FILTER CELLCHAT COMMUNICATION
# ================================================================

cellchat_AD <- filterCommunication(
  cellchat_AD,
  min.cells = 10
)

cellchat_Normal <- filterCommunication(
  cellchat_Normal,
  min.cells = 10
)

cat(
  "\nAD communications:",
  nrow(cellchat_AD@net$prob),
  "\n"
)

cat(
  "Normal communications:",
  nrow(cellchat_Normal@net$prob),
  "\n"
)

cat(
  "\nSECTION 61 COMPLETED\n"
)


# ================================================================
# SECTION 62 — CELLCHAT NETWORK VISUALIZATION
# ================================================================

cellchat_AD <- aggregateNet(cellchat_AD)
cellchat_Normal <- aggregateNet(cellchat_Normal)

dir.create(
  "figures/CellChat",
  showWarnings = FALSE,
  recursive = TRUE
)

# Cell numbers

groupSize_AD <- as.numeric(
  table(cellchat_AD@idents)
)

groupSize_Normal <- as.numeric(
  table(cellchat_Normal@idents)
)

# AD — number of interactions

pdf(
  "figures/CellChat/AD_interaction_count_network.pdf",
  width = 10,
  height = 10
)

netVisual_circle(
  cellchat_AD@net$count,
  vertex.weight = groupSize_AD,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Alzheimer's Disease - Interaction Count"
)

dev.off()

# AD — interaction strength

pdf(
  "figures/CellChat/AD_interaction_strength_network.pdf",
  width = 10,
  height = 10
)

netVisual_circle(
  cellchat_AD@net$weight,
  vertex.weight = groupSize_AD,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Alzheimer's Disease - Interaction Strength"
)

dev.off()

# Normal — number of interactions

pdf(
  "figures/CellChat/Normal_interaction_count_network.pdf",
  width = 10,
  height = 10
)

netVisual_circle(
  cellchat_Normal@net$count,
  vertex.weight = groupSize_Normal,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Normal - Interaction Count"
)

dev.off()

# Normal — interaction strength

pdf(
  "figures/CellChat/Normal_interaction_strength_network.pdf",
  width = 10,
  height = 10
)

netVisual_circle(
  cellchat_Normal@net$weight,
  vertex.weight = groupSize_Normal,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Normal - Interaction Strength"
)

dev.off()

cat(
  "\nSECTION 62 COMPLETED\n"
)

cat(
  "CellChat network figures saved in figures/CellChat/\n"
)


# ================================================================
# SECTION 63 — AD VS NORMAL COMMUNICATION COMPARISON
# ================================================================

cellchat_list <- list(
  AD = cellchat_AD,
  Normal = cellchat_Normal
)

cellchat_merged <- mergeCellChat(
  cellchat_list,
  add.names = names(cellchat_list)
)

# Compare total number of interactions

pdf(
  "figures/CellChat/AD_vs_Normal_interaction_count.pdf",
  width = 8,
  height = 6
)

compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "count"
)

dev.off()

# Compare interaction strength

pdf(
  "figures/CellChat/AD_vs_Normal_interaction_strength.pdf",
  width = 8,
  height = 6
)

compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "weight"
)

dev.off()

cat(
  "\nSECTION 63 COMPLETED\n"
)

cat(
  "AD vs Normal communication comparison saved.\n"
)
# ================================================================
# SECTION 64 — COMPUTE SIGNALING PATHWAYS
# ================================================================

cellchat_AD <- computeCommunProbPathway(
  cellchat_AD
)

cellchat_Normal <- computeCommunProbPathway(
  cellchat_Normal
)

cellchat_AD <- aggregateNet(
  cellchat_AD
)

cellchat_Normal <- aggregateNet(
  cellchat_Normal
)

cellchat_merged <- mergeCellChat(
  list(
    AD = cellchat_AD,
    Normal = cellchat_Normal
  ),
  add.names = c(
    "AD",
    "Normal"
  )
)

pathway_comparison <- rankNet(
  cellchat_merged,
  mode = "comparison",
  stacked = TRUE,
  do.stat = TRUE
)

pathway_comparison

ggsave(
  "figures/CellChat/AD_vs_Normal_signaling_pathways.png",
  pathway_comparison,
  width = 12,
  height = 8
)

cat(
  "\nSECTION 64 COMPLETED\n"
)


# ================================================================
# SECTION 65 — SAVE CELLCHAT RESULTS
# ================================================================

dir.create(
  "data/processed/CellChat",
  showWarnings = FALSE,
  recursive = TRUE
)

saveRDS(
  cellchat_AD,
  "data/processed/CellChat/GSE243292_CellChat_AD.rds"
)

saveRDS(
  cellchat_Normal,
  "data/processed/CellChat/GSE243292_CellChat_Normal.rds"
)

saveRDS(
  cellchat_merged,
  "data/processed/CellChat/GSE243292_CellChat_merged.rds"
)

cat(
  "\n============================================\n"
)

cat(
  "SECTION 65 COMPLETED\n"
)

cat(
  "CellChat AD object saved.\n"
)

cat(
  "CellChat Normal object saved.\n"
)

cat(
  "Merged CellChat object saved.\n"
)

cat(
  "============================================\n"
)



# ================================================================
# SECTION 66 — FINAL PROJECT CHECKPOINT
# ================================================================

cat(
  "\n============================================\n"
)

cat(
  "GSE243292 SINGLE-NUCLEUS RNA-SEQ ANALYSIS\n"
)

cat(
  "FINAL CHECKPOINT\n"
)

cat(
  "============================================\n\n"
)

cat(
  "Cells:",
  ncol(seurat_obj),
  "\n"
)

cat(
  "Genes:",
  nrow(seurat_obj),
  "\n"
)

cat(
  "Clusters:",
  length(
    unique(
      seurat_obj$seurat_clusters
    )
  ),
  "\n"
)

cat(
  "Cell types:",
  length(
    unique(
      na.omit(
        seurat_obj$cell_type
      )
    )
  ),
  "\n\n"
)

cat(
  "Disease groups:\n"
)

print(
  table(
    seurat_obj$Disease_Status
  )
)

cat(
  "\nMajor analyses completed:\n"
)

cat(
  "1. Quality control\n"
)

cat(
  "2. Dimensionality reduction and clustering\n"
)

cat(
  "3. Cell-type annotation\n"
)

cat(
  "4. Marker gene identification\n"
)

cat(
  "5. Disease-status analysis\n"
)

cat(
  "6. Pseudobulk differential expression\n"
)

cat(
  "7. GO functional enrichment\n"
)

cat(
  "8. OPC → Oligodendrocyte trajectory analysis\n"
)

cat(
  "9. Cell-cell communication analysis\n"
)

cat(
  "10. AD vs Normal signaling comparison\n"
)

cat(
  "\nResults directory:\n"
)

cat(
  "results/\n"
)

cat(
  "\nFigures directory:\n"
)

cat(
  "figures/\n"
)

cat(
  "\n============================================\n"
)

cat(
  "PROJECT ANALYSIS COMPLETED\n"
)

cat(
  "============================================\n"
)
