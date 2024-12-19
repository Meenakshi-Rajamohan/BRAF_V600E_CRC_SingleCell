
# Exploring Immune Infiltration in BRAF V600E Mutated Tumors of CRC using Single-Cell Transcriptomics Data

## Overview

This project investigates the immune landscape and cellular composition of BRAF V600E mutant colorectal cancer (CRC) tumors using single-cell RNA sequencing (scRNA-seq). By integrating data from two independent CRC studies, this analysis provides insights into the immunosuppressive tumor microenvironment (TME) in BRAF V600E tumors, highlighting key differences in immune cell populations compared to non-mutant tumors.

## Key Findings

- **T Regulatory Cells (Tregs)**: BRAF V600E tumors showed a marked increase in Tregs, indicating an immunosuppressive TME.
- **Cytotoxic CD8+ T Cells**: These cells were significantly depleted in BRAF V600E tumors, further exacerbating immunosuppression.
- **Myeloid Reprogramming**: Substantial depletion and functional changes were observed in macrophage subgroups, suggesting complex myeloid cell reprogramming.

## Methodology

1. **Dataset Preprocessing**:
   - scRNA-seq data from two CRC studies were quality controlled to exclude low-quality cells and doublets.
   - Batch effects were corrected using Seurat's anchor-based integration.

2. **Dimensionality Reduction**:
   - PCA, UMAP, and t-SNE were used to reduce dimensionality and visualize data.

3. **Cell Clustering and Annotation**:
   - Louvain clustering and the ProjecTILs package identified immune and myeloid cell sub-populations.
   - Marker genes were used to annotate clusters.

4. **Gene Set Enrichment Analysis (GSEA)**:
   - GSEA was employed to analyze macrophage gene signatures and other pathways related to the TME.

## Repository Contents

- `scripts/`: Modular R scripts for each step of the pipeline:
  - `01_quality_control.R`: Preprocesses raw scRNA-seq data and computes QC metrics.
  - `02_dimensionality_reduction.R`: Performs PCA, UMAP, and t-SNE for visualization.
  - `03_integration.R`: Integrates datasets and corrects batch effects.
  - `04_tcell_annotation.R`: Annotates T-cell subpopulations using ProjecTILs.
- `output/`: Contains visualization outputs.
