# GSE243292 Single-Nucleus RNA-seq Analysis of Alzheimer's Disease

## Overview

This project presents an end-to-end single-nucleus RNA-seq (snRNA-seq) analysis workflow for investigating cellular composition, cell-type-specific gene expression, disease-associated transcriptional changes, cellular trajectories, functional enrichment, and predicted cell-cell communication in Alzheimer's disease.

The analysis uses GEO dataset **GSE243292** and was primarily implemented in R using Seurat, edgeR, clusterProfiler, Monocle3, and CellChat.

## Biological Objective

The main objective is to characterize cellular and molecular changes associated with Alzheimer's disease at single-nucleus resolution.

The workflow investigates:

* Major cell populations
* Cell-type annotation
* Cluster-specific marker genes
* Alzheimer's Disease vs Normal transcriptional changes
* Cell-type-specific differential expression
* Sample-level pseudobulk differential expression
* Functional enrichment
* OPC-to-oligodendrocyte trajectory
* Cell-cell communication
* AD vs Normal signaling differences

## Dataset

**GEO accession:** GSE243292  
**Data type:** Single-nucleus RNA sequencing  
**Samples:** 15

Disease groups represented in the analyzed metadata:

|Disease group|Nuclei|
|-|-:|
|Alzheimer's Disease|61,552|
|Normal|19,892|
|Pathological Aging|41,162|

The final Seurat object contained **122,606 nuclei**, **26,102 genes**, **26 clusters**, and **8 annotated populations**.

## Cell Populations

* Astrocytes
* Endothelial cells
* Excitatory neurons
* Inhibitory neurons
* Microglia
* Oligodendrocytes
* OPCs
* Unresolved cells

Cell identities were assigned using cluster-specific marker expression and marker visualization.

## Analysis Workflow

### 1\. Quality Control

Single-nucleus data were processed using Seurat. Nuclei were retained if they had **200–4,000 detected genes** (`nFeature_RNA`) and **less than 5% mitochondrial reads** (`percent.mt`). Genes detected in fewer than 10 cells were removed at object creation (`min.cells = 10`). These thresholds should be reviewed against the QC plots in `figures/QC_before_filtering.png` and `figures/QC_after_filtering.png` for other datasets before reuse.

### 2\. Dimensionality Reduction and Clustering

PCA (30 components, top 2,000 highly variable genes) and UMAP (PCs 1–20) were used for dimensionality reduction and visualization. Graph-based (Louvain) clustering was performed at **resolution 0.5** on PCs 1–20 and identified 26 transcriptionally distinct clusters. A fixed random seed (`set.seed(42)`) is used before each stochastic step (PCA, clustering, UMAP, t-SNE) so cluster numbering and layout are reproducible between runs — this matters because cell-type labels are subsequently assigned by hardcoded cluster number (see Cell-Type Annotation).

**No batch correction or dataset integration (e.g., Harmony) was applied before clustering.** All 15 samples were merged and clustered directly on PCA space. This means clusters and the downstream disease comparison could in principle reflect sample-of-origin (technical) variation as well as biological variation; inspect `figures/UMAP_by_sample.png` alongside `figures/UMAP_by_disease.png` to check whether samples separate independently of cell type before trusting composition or DE results.

**No doublet detection or removal step (e.g., scDblFinder, DoubletFinder) was performed.** Standard nFeature/percent.mt filtering does not reliably remove doublets, which can appear as spurious "hybrid" transcriptional states.

### 3\. Cell-Type Annotation

Clusters were annotated using marker genes. An annotated UMAP was generated to visualize the major cellular populations.

### 4\. Marker Gene Identification

Cluster-specific markers were identified using Seurat `FindAllMarkers()` with `only.pos = TRUE`, `min.pct = 0.25`, and `logfc.threshold = 0.25`. Top marker tables and a marker heatmap were generated.

### 5\. Disease Status and Cell Composition

Disease status was incorporated into the metadata and visualized on UMAP. Cell-type composition was compared across Alzheimer's Disease, Normal, and Pathological Aging groups.

### 6\. Differential Expression

The primary disease comparison was **Alzheimer's Disease vs Normal**. Differential expression was evaluated within annotated cell populations.

A sample-level pseudobulk strategy was also used because individual nuclei should not be treated as independent biological replicates.

### 7\. Pseudobulk Differential Expression

Counts were aggregated by biological sample within each cell type and analyzed with edgeR.

The workflow included filtering, normalization, dispersion estimation, quasi-likelihood modeling, and Alzheimer's Disease vs Normal testing.

Significance criteria were **FDR ≤ 0.05** and **|log2FC| ≥ 0.25**.

|Cell type|Significant genes|Up in AD|Down in AD|
|-|-:|-:|-:|
|Excitatory neurons|24|4|20|
|Inhibitory neurons|10|0|10|
|Astrocytes|0|0|0|
|OPC|4|0|4|
|Microglia|3|0|3|
|Oligodendrocytes|3|0|3|
|Endothelial|0|0|0|

### 8\. Functional Enrichment

GO Biological Process enrichment was performed using clusterProfiler and org.Hs.eg.db. Significant upregulated and downregulated genes were analyzed separately.

GO enrichment was obtained for several populations, including excitatory neurons, inhibitory neurons, and microglia. Some populations had too few significant genes for enrichment under the selected thresholds.

### 9\. OPC-to-Oligodendrocyte Trajectory

Monocle3 was used to investigate an inferred transcriptional trajectory between OPCs and oligodendrocytes.

Trajectory populations:

* OPCs: 6,794
* Oligodendrocytes: 43,121

OPCs were used as the root population.

|Cell type|Median pseudotime|
|-|-:|
|OPC|0.027|
|Oligodendrocytes|2.251|

The inferred direction was therefore consistent with **OPC → Oligodendrocyte**.

Trajectory-associated genes were identified using Monocle3 `graph\_test()` and visualized across pseudotime.

### 10\. Cell-Cell Communication

CellChat was used to investigate predicted ligand-receptor communication between Astrocytes, Endothelial cells, Excitatory neurons, Inhibitory neurons, Microglia, Oligodendrocytes, and OPCs. Unresolved cells were excluded.

The workflow included signaling-gene identification, overexpressed gene and interaction identification, communication probability calculation, filtering, network aggregation, AD vs Normal comparison, and signaling pathway comparison.

CellChat results are stored under `data/processed/CellChat/`, with figures under `figures/CellChat/`.

## Software and Tools

|Tool|Purpose|
|-|-|
|R|Main analysis environment|
|Seurat|QC, clustering, UMAP and marker analysis|
|edgeR|Pseudobulk differential expression|
|clusterProfiler|Functional enrichment|
|org.Hs.eg.db|Human gene annotation|
|Monocle3|Trajectory and pseudotime analysis|
|CellChat|Cell-cell communication|
|ggplot2|Visualization|
|Harmony|Available for optional batch integration — **not used in this run**; clustering was performed on merged, non-integrated PCA space (see Quality Control / Dimensionality Reduction above)|

## Project Structure

```text
GSE243292-snRNAseq/
│
├── README.md
├── biological_interpretation.md
├── .gitignore
│
├── python/
│   └── 01_convert_h5ad.py
│
├── R/
│   └── sn_RNA_Alzhimer.R
│
├── metadata/
│   ├── sample_metadata.csv
│   ├── cell_metadata.tsv
│   ├── gene_metadata.tsv
│   └── GEO_sample_metadata/
│       ├── GSM7782916
│       ├── GSM7782917
│       ├── GSM7782918
│       ├── ...
│       └── GSM7782930
│
├── figures/
│   ├── QC/
│   ├── PCA/
│   ├── UMAP/
│   ├── annotation/
│   ├── disease/
│   ├── CellChat/
│   └── trajectory/
│
├── results/
│   ├── marker_genes/
│   ├── cell_composition/
│   ├── differential_expression/
│   ├── pseudobulk_DE/
│   ├── pathway_enrichment/
│   ├── trajectory/
│   └── CellChat/
```
**Note:** Analysis outputs are organized into dedicated subdirectories within `results/` and `figures/` for better accessibility and reproducibility. Trajectory-related results, including pseudotime analysis and significant trajectory-associated genes, are stored under `results/trajectory/`, with corresponding visualizations under `figures/trajectory/`. CellChat summary results are maintained under `results/CellChat/`, while communication network visualizations are stored under `figures/CellChat/`. Large processed Seurat and CellChat `.rds` objects are stored externally on Google Drive due to their file size. Marker gene tables, cell-type composition results, differential expression results, and pathway enrichment outputs are organized into their respective directories.
## Data

This project uses the Alzheimer's disease single-nucleus RNA-seq
dataset GSE243292.

Due to the large size of the raw, converted, and processed
single-cell datasets, these files are not stored directly in
this GitHub repository.

### Dataset

GEO accession:

GSE243292

### Data files

The following large files are stored externally:

- GSE243292_ADsnRNAseq_GEO.h5ad
- GSE243292_ADsnRNAseq_GEO.h5seurat
- Processed Seurat objects (.rds)
- Large processed expression matrices
- Large CellChat/intermediate analysis objects

### Download

Large processed files:

[Google Drive – GSE243292 processed data]([GOOGLE_DRIVE_LINK](https://drive.google.com/drive/folders/16jRuk9QApE5ozmVHT50An3NH_HkAaJHH?usp=drive_link))

After downloading, place the files according to the directory
structure used by the analysis scripts.

## Major Outputs

* QC results
* Cluster assignments
* Annotated UMAP
* Disease-status UMAP
* Marker tables and heatmaps
* Cell-type composition tables
* Pseudobulk DE tables
* Significant DE gene lists
* GO enrichment results
* OPC-to-oligodendrocyte trajectory
* Pseudotime values
* Trajectory-associated genes
* CellChat communication networks
* AD vs Normal communication comparisons
* Signaling pathway comparisons

## Reproducibility

Important processed objects are saved as RDS files:

```text
data/processed/GSE243292\_disease\_DE\_Seurat.rds
data/processed/GSE243292\_OPC\_Oligodendrocyte\_trajectory.rds
data/processed/CellChat/GSE243292\_CellChat\_AD.rds
data/processed/CellChat/GSE243292\_CellChat\_Normal.rds
data/processed/CellChat/GSE243292\_CellChat\_merged.rds
```

## Limitations

* Single-nucleus data represent nuclear rather than whole-cell transcriptomes.
* Cell-type annotation depends on marker-gene expression and may contain uncertainty.
* Monocle3 pseudotime represents an inferred transcriptional ordering rather than experimentally demonstrated lineage tracing.
* CellChat predicts ligand-receptor communication and does not experimentally validate interactions.
* Pseudobulk DE provides a more appropriate framework for biological replication than treating individual nuclei as independent samples.
* Some cell types contained too few significant genes for reliable enrichment analysis.
* **The Alzheimer's Disease vs Normal comparison is based on 8 AD samples vs only 2 Normal samples.** This is a small biological-replicate count, particularly for the Normal group, and limits statistical power to detect anything but relatively large, consistent effects. Non-significant results for a given cell type should not be read as evidence of "no change."
* **The pseudobulk DE model (`~ Disease_Status`) does not include covariates** such as age, sex, or postmortem interval (PMI), which are standard potential confounders in postmortem brain transcriptomics and were not available/modeled in this analysis. Some disease-associated genes could partly reflect unmodeled donor-level differences rather than disease status alone.
* **No batch correction or dataset integration (e.g., Harmony) was applied before clustering** — all 15 samples were clustered together in a single merged PCA space. Sample-of-origin technical variation could in principle influence cluster boundaries and, downstream, cell-type composition and pseudobulk DE results.
* **No doublet detection/removal step was performed.** Standard `nFeature_RNA`/`percent.mt` filtering does not reliably identify doublets, which could contribute to some of the smaller or more ambiguous clusters (e.g., "Unresolved").
* Additional donor-aware modeling and experimental validation could strengthen the conclusions.

## Conclusion

This project implements a comprehensive single-nucleus RNA-seq workflow for Alzheimer's disease using GSE243292.

The workflow integrates cellular clustering and annotation, disease-associated transcriptional analysis, sample-level pseudobulk differential expression, functional enrichment, trajectory inference, and predicted cell-cell communication.

Together, these analyses provide an integrated view of cellular composition, molecular changes, oligodendrocyte-lineage states, and intercellular signaling associated with Alzheimer's disease.

## Acknowledgement

This project uses publicly available data from the Gene Expression Omnibus (GEO), accession **GSE243292**.
