# Beginner Spatial Transcriptomics Workflow in R

This repository contains a beginner-friendly R workflow for analyzing 10x Genomics Visium spatial transcriptomics data using Seurat. The workflow is designed around four lacrimal gland tissue sections: two control samples and two disease-model samples. It focuses on quality control, normalization, integration, clustering, spatial visualization, and marker-based cluster interpretation.

This script does not perform disease-versus-control differential expression analysis. Instead, it uses marker genes and tissue localization to support biological interpretation of spatial clusters.

---

## Project Goal

The goal of this project is to provide a clear, modular, and reproducible spatial transcriptomics workflow for beginners who want to analyze Visium data while preserving spatial context.

The workflow helps answer questions such as:

- Are the Visium sections usable based on QC metrics?
- How do the four samples integrate into a shared transcriptional space?
- Which clusters are detected after integration?
- Where are those clusters located in the tissue?
- Which marker genes help interpret the biological identity of each cluster?

---

## Original Study and Data Source

This workflow is based on the publicly available lacrimal gland spatial transcriptomics dataset from:

**Mauduit, O., Delcroix, V., Umazume, T., de Paiva, C. S., Dartt, D. A., & Makarenkova, H. P. (2022). _Spatial transcriptomics of the lacrimal gland features macrophage activity and epithelium metabolism as key alterations during chronic inflammation_. Frontiers in Immunology, 13, 1011125. https://doi.org/10.3389/fimmu.2022.1011125**

The original authors deposited the study data in the **NCBI Gene Expression Omnibus (GEO)**:

- **Bulk RNA-seq dataset:** GSE210332
- **Visium spatial gene expression dataset:** GSE210380

The original article also provides supplemental files, including a ready-to-use Seurat object and analysis code. This repository is a beginner-oriented re-analysis and teaching workflow based on those publicly available resources. The original data, biological study design, and primary scientific findings belong to the original authors.

Please cite the original Frontiers in Immunology paper if you use this workflow or the associated dataset.

---

## Dataset Design

The workflow is organized around four Visium sections.

| Project Label | Original Sample ID | Condition | Approximate Reads |
|---|---|---:|---:|
| Control_1 | BALB/c_1 | Control | 49 million |
| Control_2 | BALB/c_2 | Control | 380 million |
| Case_1 | NOD.H-2b_1 | Case | 288 million |
| Case_2 | NOD.H-2b_2 | Case | 101 million |

`Control_2` and `Case_1` are used as representative sections for spatial cluster overlays and spatial marker-expression plots because they are the deeper-sequenced control and case samples.

---

## Main R Packages

The workflow uses the following R packages:

```r
Seurat
ggplot2
patchwork
dplyr
tidyr
future
cowplot
png
```

The script checks whether these packages are installed before running. If any package is missing, the script stops and prints the package names that need to be installed.

---

## Expected Input Folder Structure

Before running the script, edit the `root_dir` path so it points to the folder containing your four sample folders.

Example:

```text
spatial transcriptomics/
├── BALBC_1/
│   └── processed data/
│       ├── filtered_feature_bc_matrix.h5
│       └── spatial/
├── BALBC_2/
│   └── processed data/
│       ├── filtered_feature_bc_matrix.h5
│       └── spatial/
├── NOD_1/
│   └── processed data/
│       ├── filtered_feature_bc_matrix.h5
│       └── spatial/
└── NOD_2/
    └── processed data/
        ├── filtered_feature_bc_matrix.h5
        └── spatial/
```

The script can detect Space Ranger outputs either directly inside each sample folder or inside a subfolder named `processed data`.

---

## How to Run the Workflow

1. Open the R script in RStudio.
2. Update the project path:

```r
root_dir <- "C:/Users/yourname/path/to/spatial transcriptomics"
```

3. Confirm that each sample folder contains:

```text
filtered_feature_bc_matrix.h5
spatial/
```

4. Run the script from top to bottom.

All output files will be saved automatically in:

```text
beginner_spatial_outputs_no_DEG/
```

---

## Workflow Modules

### 1. Project Setup

Defines required packages, sample metadata, memory settings, analysis parameters, and output folders.

### 2. Figure Formatting

Standardizes plot fonts, figure sizes, and 300 dpi export settings so the figures can be inserted into reports, manuscripts, or Word documents.

### 3. Helper Functions

Includes reusable functions for:

- Finding Space Ranger output folders
- Loading Visium samples
- Reading H&E tissue images
- Saving plots
- Saving tables

### 4. Loading Visium Samples

Loads each sample into a separate Seurat object using `Load10X_Spatial()` and adds sample-level metadata such as condition and original sample ID.

### 5. Quality Control

Calculates and exports spot-level QC summaries, including:

- Number of spots
- Median UMI count
- Median detected features
- Median mitochondrial percentage

QC violin plots are generated for each sample.

### 6. SCTransform Normalization

Normalizes each sample independently using `SCTransform()`.

### 7. Integration

Integrates the four normalized samples into one shared Seurat object using SCT-based integration.

### 8. PCA, UMAP, and Clustering

Performs dimensionality reduction and clustering using:

- PCA
- UMAP
- Shared nearest-neighbor clustering

Cluster labels are numeric at this stage. Biological interpretation is performed later using marker genes and spatial context.

### 9. Marker Gene Analysis

Identifies positive marker genes for each cluster using `FindAllMarkers()`.

This step is used for cluster interpretation only. It is not a case-control differential expression analysis.

### 10. UMAP Visualization

Generates:

- Integrated UMAP showing all clusters
- Split UMAP showing cluster distribution across individual samples

### 11. Spatial Cluster Overlay

Projects cluster identities back onto the tissue sections for the representative control and case samples.

### 12. Marker Expression Visualization

Plots selected marker genes on both UMAP space and tissue space.

Marker genes used in this workflow include:

| Marker | Suggested Signal |
|---|---|
| Scgb2b17 | Acinar-associated epithelial signal |
| Ltf | Ductal-associated epithelial signal |
| Jchain | Plasma-cell-rich signal |
| Cd79a | B-cell-associated signal |
| Cd3g | T-cell-associated signal |
| Gsn | Macrophage or stromal-associated signal |

---

## Main Output Files

### Seurat Objects

| File | Description |
|---|---|
| `01_raw_visium_objects.rds` | Raw Seurat objects after loading Visium data |
| `02_sct_normalized_objects.rds` | SCTransform-normalized sample objects |
| `03_integrated_object_raw.rds` | Integrated object before clustering |
| `04_integrated_object_clustered.rds` | Integrated object after PCA, UMAP, and clustering |
| `final_integrated_spatial_object_no_DEG.rds` | Final processed Seurat object |

### Tables

| File | Description |
|---|---|
| `sample_metadata_mapping.csv` | Sample labels, conditions, and sequencing-depth metadata |
| `QC_summary_by_sample.csv` | QC summary for each Visium section |
| `cluster_markers_all.csv` | All positive cluster marker genes |
| `cluster_markers_top10.csv` | Top 10 marker genes per cluster |
| `cluster_interpretation_marker_guide.csv` | Marker guide for biological interpretation |

### Figures

| File | Description |
|---|---|
| `QC_violin_plots_all_samples.pdf` | QC violin plots for all samples |
| `QC_violin_plot_<sample>.png` | Sample-specific QC violin plots |
| `plot_4A_integrated_umap.png` | Integrated UMAP |
| `plot_4A_split_umap.png` | UMAP split by sample |
| `plot_4B_spatial_clusters.png` | Spatial cluster overlays on representative tissue sections |
| `plot_4C_feature_plots.png` | Marker FeaturePlots on UMAP |
| `plot_spatial_marker_expression_control.png` | Spatial marker expression in representative control |
| `plot_spatial_marker_expression_case.png` | Spatial marker expression in representative case |

---

## Important Notes

- Update `root_dir` before running the script.
- The workflow assumes Space Ranger output files are already generated.
- The script uses `future::plan("sequential")` for stability on standard laptops and desktops.
- The memory limit is set to 8 GB using `future.globals.maxSize`.
- Cluster numbers may change if parameters such as resolution, integration features, or PCA dimensions are modified.
- Marker-based cluster interpretation should be treated as biological guidance, not as fixed annotation.
- Low-count spots in Visium data are not always poor-quality spots. They may reflect tissue boundaries, immune-rich regions, or regions with lower RNA content.
- Large raw data files and Seurat `.rds` objects should usually not be uploaded directly to GitHub. Use GEO, institutional storage, cloud storage, Zenodo, Figshare, or Git LFS for large files.

---

## GitHub Description

Beginner-friendly Seurat workflow for 10x Visium spatial transcriptomics analysis of lacrimal gland tissue, including QC, SCT integration, clustering, UMAP visualization, spatial overlays, and marker-based cluster interpretation.

---

## Citation


