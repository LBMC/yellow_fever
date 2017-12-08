---
title: "Quality Control Analysis"
author: Laurent Modolo [laurent.modolo@ens-lyon.fr](mailto:laurent.modolo@ens-lyon.fr)
date: 2017/11/22
output:
  pdf_document:
    toc: true
    toc_depth: 1
    number_section: true
    highlight: tango
    latex_engine: xelatex
---

Every step of this analysis is wrapped inside the package `scRNAtools`.
To run this analysis in R, you can run the script `src/2_cell_type.R` after installing
this package.

# Manual annotation of the cell-type

T-cells are classically classified as **CM**, **TSCM**, **EM** or **TEMRA**.
Throught this document we define **MEM** as **CM** or **TSCM** looking cells
and **EFF** as **EM** or **TEMRA** looking cells. **CM** and **TSCM** beeing
more of the memory type and **EM** and **TEMRA** beeing more of the effector
type.

We start with the following set of manually annotated cells (based on *ccr7*
and *cd45ra*):

|CM | TSCM | EM | TEMRA
----|------|----|------
|61 | 40   | 43 | 58

Which gives us the following table for *phenotype_surface_marker*:

|MEM | EFF
|----|---
|101 | 101

we load the data with the following commands:
```R
load("results/QC/CB_counts_QC.Rdata")
scd$setfeature("surface_cell_type", phenotype_surface_marker)
b_cells <- scd$getfeature("QC_good") %in% T
```

# PLS classification based the manual annotation

In the following we want to extend the manual annotation of a subset of
cells in the male donor at day 15, 136 and 593 to all cells for those 3 time
points.

Instead of on relying on the *ccr7* and *cd45ra* surface markers used to
classify the 4 T-cells types. We use the following surface markers and genes
to make our classification based on the literature.

```R
surface_marker <- c(
"ccr7", "il7ra"
)
genes_marker <- c(
"GNLY", "GZMH", "CCL4", "KLRD1" "GZMB", "ZEB2", "LTB",  "TCF7", "CCR7"
"GZMK", "SELL"
)
```

We then fit a sparse logistic PLS model on the annotated subset of cells using
these features with the following command:

```R
surface_cell_type_classification <- classification(
  scd = scd$select(b_cells = b_cells),
  feature = "surface_cell_type",
  features = surface_marker,
  genes = genes_marker,
  ncores = 10,
  algo = "spls_stab",
  output_file = "results/cell_type/surface_cell_types"
)
save(
  surface_cell_type_classification,
  file = "results/cell_type/surface_cell_types_all_smplscv.Rdata"
)
cell_type_groups <- rep(NA, scd$getncells)
cell_type_groups[b_cells] <- cell_type_classification$groups
scd$setfeature("surface_cell_type", cell_type_groups)
save(scd, file = "results/cell_type/CB_counts_QC_surface_cell_type.Rdata")
```

We also export this new feature to the raw counts:

```R
scd_norm <- scd
load("results/QC/cells_counts_QC.Rdata")
scd <- scdata$new(
  infos = scd_norm$getfeatures,
  counts = scd$getcounts
)
save(scd, file = "results/cell_type/cells_counts_QC_surface_cell_type.Rdata")
```

\begin{center}
  \begin{figure}
    \includegraphics[width=0.5\textwidth]{../results/cell_type/pca/pca_counts_QC_surface_cell_type.pdf}\includegraphics[width=0.5\textwidth]{../results/cell_type/pcmf/pcmf_counts_QC_surface_cell_type.pdf}
    \caption{PCA and pCMF plot for cells classification on all genes and all day)}
  \end{figure}
\end{center}

\begin{center}
  \begin{figure}
    \includegraphics[width=0.3\textwidth]{../results/cell_type/pca/pca_counts_QC_surface_cell_type_D15.pdf} \includegraphics[width=0.3\textwidth]{../results/cell_type/pca/pca_counts_QC_surface_cell_type_D136.pdf} \includegraphics[width=0.3\textwidth]{../results/cell_type/pca/pca_counts_QC_surface_cell_type_D593.pdf}
    \caption{PCA plot for cells classification on all genes)}
  \end{figure}
\end{center}

\begin{center}
  \begin{figure}
    \includegraphics[width=0.3\textwidth]{../results/cell_type/pcmf/pcmf_counts_QC_surface_cell_type_D15.pdf} \includegraphics[width=0.3\textwidth]{../results/cell_type/pcmf/pcmf_counts_QC_surface_cell_type_D136.pdf} \includegraphics[width=0.3\textwidth]{../results/cell_type/pcmf/pcmf_counts_QC_surface_cell_type_D593.pdf}
    \caption{pCMF plot for cells classification on all genes)}
  \end{figure}
\end{center}

# Differential expression analysis

To improve our classification of the different cell-type
