# MorPhiC Bulk RNA-seq Explorer

A Shiny application for browsing, comparing, and analysing differential gene expression across MorPhiC knock-out assays.

## Pages

| Page | Description |
|------|-------------|
| **Home** | Landing page with quick-access cards to each section |
| **Perturbed Genes Overview** | Study-level summaries, functional annotations (Reactome, GO, PANTHER), and disease/phenotype associations (IMPC, HPO, OMIM, Orphanet) |
| **Assay Comparison** | Multi-assay side-by-side DEG matrix, UpSet intersection plots, downloadable gene lists |
| **Assay Analysis** | Single-assay deep dive with pathway overlaps, GO terms, protein classes, and phenotype/disease browsers |

## Data

All data is served from a local **DuckDB** database (`morphic_bulkRNA3.duckdb`) containing:

- **170+ assay tables** — per-gene DEG status and log2 fold-change
- **study_info** — assay metadata (gene, model system, KO strategy, cell line, DPC, etc.)
- **Reference tables** — Reactome, GO BP/MF, PANTHER class/family gene sets
- **Annotation tables** — IMPC, HPO, OMIM, Orphanet phenotype/disease mappings

## Modules

```
modules/
  helpers.R          # Shared constants, colour palettes, utility functions
  home.R             # Landing page
  studies_overview.R # Perturbed Genes Overview
  compare_studies.R  # Assay Comparison
  deg_annotations.R  # Assay Analysis (annotations, phenotypes, diseases)
```

## Requirements

R 4.x with the following packages:

```
shiny, bslib, bsicons, DT, plotly, duckdb, DBI, dplyr, UpSetR, shinyWidgets
```

Install all dependencies:

```r
source("dependencies.R")
```

## Running

```r
shiny::runApp(".")
```

## Colour conventions

| DEG direction | Colour |
|---------------|--------|
| Up-regulated | `#E63946` (red) |
| Down-regulated | `#457B9D` (blue) |
| All / combined | `#6A4C93` (purple) |
