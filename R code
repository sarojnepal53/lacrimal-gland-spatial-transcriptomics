############################################################
# Beginner Spatial Transcriptomics Guide in R
# Modular Seurat workflow for 10x Visium lacrimal gland data
#
# Focus:
# - Four Visium sections
# - Two controls and two disease-model samples
# - No disease-versus-control DEG module
# - Cluster interpretation is based on marker genes and spatial context
#
# Original case-study sample mapping:
# - Control_1 = BALB/c_1
# - Control_2 = BALB/c_2
# - Case_1    = NOD.H-2b_1
# - Case_2    = NOD.H-2b_2
#
# Sequencing-depth context reported in the paper:
# - BALB/c_1:     about 49 million reads
# - BALB/c_2:     about 380 million reads
# - NOD.H-2b_1:   about 288 million reads
# - NOD.H-2b_2:   about 101 million reads
#
# Control_2 and Case_1 were the deeper-sequenced sections,
# so this script uses them as representative sections for
# spatial overlay and spatial marker-expression plots.
############################################################


############################################################
# Module 1: Project setup, sample metadata, and memory handling
############################################################

# This module sets up the project folder, required packages,
# sample labels, biological conditions, and memory-safe settings.

required_pkgs <- c(
  "Seurat",
  "ggplot2",
  "patchwork",
  "dplyr",
  "tidyr",
  "future",
  "cowplot",
  "png"
)

missing_required <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_required) > 0) {
  stop(
    "Please install these packages first: ",
    paste(missing_required, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(future)
  library(cowplot)
  library(png)
})

set.seed(1234)

# Running sequentially is often more stable for beginners,
# especially on laptops or desktops with limited RAM.
future::plan("sequential")

# Increase maximum object size for operations such as SCT integration.
# Adjust this value if your computer has less or more RAM.
options(future.globals.maxSize = 8 * 1024^3)

# EDIT THIS PATH.
# This folder should contain the four sample folders.
# Example structure:
# C:/Users/sawro/OneDrive/Desktop/spatial transcriptomics/
#   BALBC_1/processed data/
#   BALBC_2/processed data/
#   NOD_1/processed data/
#   NOD_2/processed data/
root_dir <- "C:/Users/sawro/OneDrive/Desktop/spatial transcriptomics"

# The folder_name column uses your local folder names.
# The original_sample_id column records the paper-style labels.
# The sample_id column gives beginner-friendly labels used in plots.
sample_info <- data.frame(
  folder_name = c("BALBC_1", "BALBC_2", "NOD_1", "NOD_2"),
  original_sample_id = c("BALB/c_1", "BALB/c_2", "NOD.H-2b_1", "NOD.H-2b_2"),
  sample_id = c("Control_1", "Control_2", "Case_1", "Case_2"),
  condition = c("Control", "Control", "Case", "Case"),
  approximate_reads_million = c(49, 380, 288, 101),
  mean_reads_per_spot_reported = c(43487, 343286, 188827, 67466),
  median_genes_per_spot_reported = c(874, 2227, 1848, 1051),
  stringsAsFactors = FALSE
)

# Output folder for plots and tables.
out_dir <- file.path(root_dir, "beginner_spatial_outputs_no_DEG")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Save sample metadata for record keeping.
write.csv(
  sample_info,
  file = file.path(out_dir, "sample_metadata_mapping.csv"),
  row.names = FALSE
)

# Analysis parameters.
integration_nfeatures <- 2000
dims_use <- 1:30
cluster_resolution <- 0.42

# Representative sections used for tissue overlay and marker plots.
# These are the deeper-sequenced control and case sections.
compare_control <- "Control_2"
compare_case <- "Case_1"

# Marker genes used for biological interpretation.
# These markers follow the lacrimal gland case-study logic.
feature_panel <- c(
  "Scgb2b17",  # acinar-associated
  "Jchain",   # plasma-cell-rich spots
  "Cd79a",    # B-cell-associated
  "Cd3g",     # T-cell-associated
  "Gsn",      # macrophage/stromal-associated
  "Ltf"       # ductal-associated
)


############################################################
# Module 1B: Manuscript font size and export settings
############################################################

# IMPORTANT:
# Use "sans" instead of "Arial" for the base pdf() device.
# On most Windows systems, "sans" maps visually to Arial-like text.
# This avoids the error: unknown family 'Arial'.
# If you specifically need Arial in a PDF, use a Cairo-based device or
# configure fonts separately; do not pass family = "Arial" to base pdf().
FIG_FONT_FAMILY <- "sans"

# Font sizes suitable for manuscript figures after inserting into Word.
FIG_BASE_SIZE <- 14
FIG_TITLE_SIZE <- 16
FIG_AXIS_TITLE_SIZE <- 14
FIG_AXIS_TEXT_SIZE <- 12
FIG_LEGEND_TITLE_SIZE <- 12
FIG_LEGEND_TEXT_SIZE <- 11
FIG_STRIP_TEXT_SIZE <- 13
SEURAT_LABEL_SIZE <- 5
HE_LABEL_SIZE <- 15

# Export resolution. 300 dpi is standard for manuscript raster figures.
FIG_DPI <- 300

# Keep the same plot layout and panel dimensions as your original script.
# Only the font size and export dpi are standardized.
FIG_SINGLE_WIDTH <- 7
FIG_SINGLE_HEIGHT <- 5
FIG_QC_WIDTH <- 14
FIG_QC_HEIGHT <- 4
FIG_SPLIT_WIDTH <- 10
FIG_SPLIT_HEIGHT <- 8
FIG_SQUARE_WIDTH <- 10
FIG_SQUARE_HEIGHT <- 10
FIG_TALL_WIDTH <- 10
FIG_TALL_HEIGHT <- 12

manuscript_theme <- theme(
  text = element_text(family = FIG_FONT_FAMILY, size = FIG_BASE_SIZE),
  plot.title = element_text(
    family = FIG_FONT_FAMILY,
    size = FIG_TITLE_SIZE,
    face = "bold",
    hjust = 0.5
  ),
  axis.title = element_text(
    family = FIG_FONT_FAMILY,
    size = FIG_AXIS_TITLE_SIZE
  ),
  axis.text = element_text(
    family = FIG_FONT_FAMILY,
    size = FIG_AXIS_TEXT_SIZE
  ),
  legend.title = element_text(
    family = FIG_FONT_FAMILY,
    size = FIG_LEGEND_TITLE_SIZE
  ),
  legend.text = element_text(
    family = FIG_FONT_FAMILY,
    size = FIG_LEGEND_TEXT_SIZE
  ),
  strip.text = element_text(
    family = FIG_FONT_FAMILY,
    size = FIG_STRIP_TEXT_SIZE,
    face = "bold"
  )
)


############################################################
# Module 2: Helper functions for loading data, reading tissue
# images, and saving plots
############################################################

# This module creates reusable functions so the main workflow
# stays clean and readable.

# This function checks whether the Space Ranger output is directly
# inside the sample folder or inside a subfolder called "processed data".
find_spaceranger_dir <- function(root_dir, folder_name) {
  direct_dir <- file.path(root_dir, folder_name)
  processed_dir <- file.path(root_dir, folder_name, "processed data")
  
  direct_h5 <- file.path(direct_dir, "filtered_feature_bc_matrix.h5")
  processed_h5 <- file.path(processed_dir, "filtered_feature_bc_matrix.h5")
  
  if (file.exists(direct_h5)) {
    return(direct_dir)
  }
  
  if (file.exists(processed_h5)) {
    return(processed_dir)
  }
  
  stop(
    "Could not find filtered_feature_bc_matrix.h5 for folder: ",
    folder_name,
    "\nChecked:\n",
    direct_h5,
    "\n",
    processed_h5
  )
}

# Load one Visium sample into Seurat.
load_visium_sample <- function(sample_dir, sample_label, condition_label, original_label) {
  h5_file <- file.path(sample_dir, "filtered_feature_bc_matrix.h5")
  spatial_dir <- file.path(sample_dir, "spatial")
  
  if (!file.exists(h5_file)) {
    stop("Missing file: ", h5_file)
  }
  
  if (!dir.exists(spatial_dir)) {
    stop("Missing spatial folder: ", spatial_dir)
  }
  
  obj <- Load10X_Spatial(
    data.dir = sample_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = sample_label
  )
  
  obj$sample_id <- sample_label
  obj$condition <- condition_label
  obj$original_sample_id <- original_label
  obj$orig.ident <- sample_label
  
  # Mouse mitochondrial genes usually start with "mt-".
  # Human mitochondrial genes usually start with "MT-".
  mt_mouse <- grep("^mt-", rownames(obj), value = TRUE)
  mt_human <- grep("^MT-", rownames(obj), value = TRUE)
  
  if (length(mt_mouse) > 0) {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  } else if (length(mt_human) > 0) {
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  } else {
    obj$percent.mt <- 0
  }
  
  return(obj)
}

# Read and plot the low-resolution H&E image from a Space Ranger folder.
plot_he_image <- function(sample_dir, sample_label) {
  img_path <- file.path(sample_dir, "spatial", "tissue_lowres_image.png")
  
  if (!file.exists(img_path)) {
    return(
      ggplot() +
        theme_void() +
        ggtitle(paste(sample_label, "H&E image not found"))
    )
  }
  
  img <- png::readPNG(img_path)
  
  cowplot::ggdraw() +
    cowplot::draw_image(img) +
    cowplot::draw_label(
      sample_label,
      x = 0.5,
      y = 0.98,
      fontface = "bold",
      fontfamily = FIG_FONT_FAMILY,
      size = HE_LABEL_SIZE
    ) +
    theme_void()
}

# Save plots consistently at 300 dpi.
save_plot <- function(filename, plot_obj, width = 8, height = 6, dpi = FIG_DPI) {
  ggsave(
    filename = file.path(out_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
}

# Save a data frame consistently.
save_table <- function(filename, table_obj) {
  write.csv(
    table_obj,
    file = file.path(out_dir, filename),
    row.names = FALSE
  )
}


############################################################
# Module 3: Loading all four Visium sections into Seurat
############################################################

# This module loads each Space Ranger output folder into Seurat.
# At this stage, each section is still an independent Seurat object.

obj_list <- list()

for (i in seq_len(nrow(sample_info))) {
  current_dir <- find_spaceranger_dir(
    root_dir = root_dir,
    folder_name = sample_info$folder_name[i]
  )
  
  message("Loading sample: ", sample_info$sample_id[i])
  message("Using directory: ", current_dir)
  
  obj_list[[sample_info$sample_id[i]]] <- load_visium_sample(
    sample_dir = current_dir,
    sample_label = sample_info$sample_id[i],
    condition_label = sample_info$condition[i],
    original_label = sample_info$original_sample_id[i]
  )
}

saveRDS(
  obj_list,
  file = file.path(out_dir, "01_raw_visium_objects.rds")
)


############################################################
# Module 4: Quality control with violin plots
############################################################

# This module checks spot-level quality metrics.
# For Visium, low count or low feature spots are not always bad.
# They may reflect tissue boundaries, immune-rich regions, or
# lower RNA content compared with secretory epithelium.

qc_summary <- lapply(names(obj_list), function(nm) {
  obj <- obj_list[[nm]]
  
  data.frame(
    sample_id = nm,
    condition = unique(obj$condition),
    original_sample_id = unique(obj$original_sample_id),
    n_spots = ncol(obj),
    median_nCount_Spatial = median(obj$nCount_Spatial, na.rm = TRUE),
    median_nFeature_Spatial = median(obj$nFeature_Spatial, na.rm = TRUE),
    median_percent_mt = median(obj$percent.mt, na.rm = TRUE)
  )
}) %>%
  bind_rows()

save_table("QC_summary_by_sample.csv", qc_summary)

# Keep the same QC violin plots, but increase font size.
# Do not add family = FIG_FONT_FAMILY to base pdf() if FIG_FONT_FAMILY is Arial.
# Base pdf() may not know Arial and can produce: unknown family 'Arial'.
pdf(
  file = file.path(out_dir, "QC_violin_plots_all_samples.pdf"),
  width = FIG_QC_WIDTH,
  height = FIG_QC_HEIGHT
)

for (nm in names(obj_list)) {
  p_qc <- (
    VlnPlot(
      obj_list[[nm]],
      features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
      pt.size = 0.1,
      ncol = 3
    ) +
      plot_annotation(title = paste("QC:", nm))
  ) &
    manuscript_theme

  print(p_qc)

  # Also save a 300 dpi PNG copy for direct insertion into a manuscript.
  save_plot(
    paste0("QC_violin_plot_", nm, ".png"),
    p_qc,
    width = FIG_QC_WIDTH,
    height = FIG_QC_HEIGHT
  )
}

dev.off()


############################################################
# Module 5: Normalization of each section with SCTransform
############################################################

# This module normalizes each section independently.
# SCTransform reduces technical variation related to count depth
# while preserving meaningful biological structure.

obj_list_sct <- lapply(obj_list, function(x) {
  DefaultAssay(x) <- "Spatial"
  
  SCTransform(
    x,
    assay = "Spatial",
    verbose = FALSE,
    return.only.var.genes = FALSE
  )
})

saveRDS(
  obj_list_sct,
  file = file.path(out_dir, "02_sct_normalized_objects.rds")
)

gc()


############################################################
# Module 6: Integration across multiple sections
############################################################

# This module integrates all sections into one shared object.
# Integration helps reduce sample-specific technical differences
# so control and case sections can be visualized in the same
# molecular space.

integration_features <- SelectIntegrationFeatures(
  object.list = obj_list_sct,
  nfeatures = integration_nfeatures
)

obj_list_sct <- PrepSCTIntegration(
  object.list = obj_list_sct,
  anchor.features = integration_features,
  verbose = FALSE
)

anchors <- FindIntegrationAnchors(
  object.list = obj_list_sct,
  normalization.method = "SCT",
  anchor.features = integration_features,
  verbose = FALSE
)

combined <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  verbose = FALSE
)

saveRDS(
  combined,
  file = file.path(out_dir, "03_integrated_object_raw.rds")
)

rm(anchors)
gc()


############################################################
# Module 7: PCA, UMAP, and clustering
############################################################

# This module performs dimensionality reduction and clustering.
# The cluster labels generated here are numeric only.
# Biological names are assigned later using marker genes
# and spatial context.

DefaultAssay(combined) <- "integrated"

combined <- RunPCA(
  combined,
  npcs = max(dims_use),
  verbose = FALSE
)

combined <- RunUMAP(
  combined,
  reduction = "pca",
  dims = dims_use,
  verbose = FALSE
)

combined <- FindNeighbors(
  combined,
  reduction = "pca",
  dims = dims_use
)

combined <- FindClusters(
  combined,
  resolution = cluster_resolution,
  verbose = FALSE
)

Idents(combined) <- "seurat_clusters"

cluster_ids <- levels(Idents(combined))

cluster_colors <- setNames(
  grDevices::hcl.colors(length(cluster_ids), palette = "Dynamic"),
  cluster_ids
)

saveRDS(
  combined,
  file = file.path(out_dir, "04_integrated_object_clustered.rds")
)


############################################################
# Module 8: Marker genes and biological interpretation of clusters
############################################################

# This module identifies marker genes that help interpret clusters.
# This is not a disease-versus-control DEG module.
# The goal here is only to ask which genes help define each cluster.

DefaultAssay(combined) <- "SCT"

combined <- PrepSCTFindMarkers(
  combined,
  assay = "SCT",
  verbose = TRUE
)

Idents(combined) <- "seurat_clusters"

all_markers <- FindAllMarkers(
  combined,
  assay = "SCT",
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.10
)

top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()

save_table("cluster_markers_all.csv", all_markers)
save_table("cluster_markers_top10.csv", top_markers)

# A simple interpretation key based on the case-study paper.
# Your actual cluster numbers may differ, so use this as guidance,
# not as a fixed rule.
cluster_interpretation_guide <- data.frame(
  marker_gene = c("Scgb2b17", "Ltf", "Jchain", "Cd79a", "Cd3g", "Gsn"),
  suggested_signal = c(
    "Acinar-associated epithelial signal",
    "Ductal-associated epithelial signal",
    "Plasma-cell-rich signal",
    "B-cell-associated signal",
    "T-cell-associated signal",
    "Macrophage or stromal-associated signal"
  ),
  stringsAsFactors = FALSE
)

save_table("cluster_interpretation_marker_guide.csv", cluster_interpretation_guide)


############################################################
# Module 9: Integrated UMAP and split UMAP
############################################################

# This module creates UMAP plots.
# The integrated UMAP shows the shared cluster structure.
# The split UMAP shows how each sample contributes to that structure.

p_umap_integrated <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "seurat_clusters",
  cols = cluster_colors,
  raster = FALSE,
  label = TRUE,
  label.size = SEURAT_LABEL_SIZE
) +
  theme_classic(base_size = FIG_BASE_SIZE, base_family = FIG_FONT_FAMILY) +
  manuscript_theme +
  ggtitle("Integrated UMAP")

save_plot(
  "plot_4A_integrated_umap.png",
  p_umap_integrated,
  width = FIG_SINGLE_WIDTH,
  height = FIG_SINGLE_HEIGHT
)

p_umap_split <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "seurat_clusters",
  split.by = "sample_id",
  cols = cluster_colors,
  ncol = 2,
  raster = FALSE
) &
  theme_classic(base_size = FIG_BASE_SIZE, base_family = FIG_FONT_FAMILY) &
  manuscript_theme

save_plot(
  "plot_4A_split_umap.png",
  p_umap_split,
  width = FIG_SPLIT_WIDTH,
  height = FIG_SPLIT_HEIGHT
)


############################################################
# Module 10: Spatial overlays on top of the tissue image
############################################################

# This module projects cluster identity back onto the tissue.
# Control_2 and Case_1 are used here because they correspond
# to the deeper-sequenced representative sections.

control_obj <- subset(
  combined,
  subset = sample_id == compare_control
)

case_obj <- subset(
  combined,
  subset = sample_id == compare_case
)

p_control_overlay <- SpatialDimPlot(
  control_obj,
  group.by = "seurat_clusters",
  cols = cluster_colors,
  label = FALSE,
  crop = FALSE,
  pt.size.factor = 1.6
) +
  ggtitle(compare_control) +
  manuscript_theme

p_case_overlay <- SpatialDimPlot(
  case_obj,
  group.by = "seurat_clusters",
  cols = cluster_colors,
  label = FALSE,
  crop = FALSE,
  pt.size.factor = 1.6
) +
  ggtitle(compare_case) +
  manuscript_theme

control_dir <- find_spaceranger_dir(
  root_dir = root_dir,
  folder_name = sample_info$folder_name[sample_info$sample_id == compare_control]
)

case_dir <- find_spaceranger_dir(
  root_dir = root_dir,
  folder_name = sample_info$folder_name[sample_info$sample_id == compare_case]
)

p_control_he <- plot_he_image(
  sample_dir = control_dir,
  sample_label = compare_control
)

p_case_he <- plot_he_image(
  sample_dir = case_dir,
  sample_label = compare_case
)

p_spatial_pair <- ((p_control_overlay | p_case_overlay) /
  (p_control_he | p_case_he)) &
  manuscript_theme

save_plot(
  "plot_4B_spatial_clusters.png",
  p_spatial_pair,
  width = FIG_SQUARE_WIDTH,
  height = FIG_SQUARE_HEIGHT
)


############################################################
# Module 11: Marker FeaturePlots on the integrated UMAP
############################################################

# This module plots selected marker genes on the UMAP.
# FeaturePlot helps connect marker genes to transcriptional clusters.

DefaultAssay(combined) <- "SCT"

feature_panel_present <- intersect(
  feature_panel,
  rownames(combined)
)

if (length(feature_panel_present) > 0) {
  p_features <- FeaturePlot(
    combined,
    reduction = "umap",
    features = feature_panel_present,
    ncol = 2,
    order = TRUE,
    raster = FALSE
  ) &
    manuscript_theme

  save_plot(
    "plot_4C_feature_plots.png",
    p_features,
    width = FIG_TALL_WIDTH,
    height = FIG_TALL_HEIGHT
  )
} else {
  message("None of the selected marker genes were found in the dataset.")
}


############################################################
# Module 11B: Spatial marker expression on representative
# control and case sections
############################################################

# This module plots selected marker genes directly on tissue sections.
# SpatialFeaturePlot helps show where marker signals occur in tissue space.

DefaultAssay(control_obj) <- "SCT"
DefaultAssay(case_obj) <- "SCT"

marker_panel_present <- intersect(
  feature_panel,
  rownames(combined)
)

if (length(marker_panel_present) > 0) {
  control_marker_spatial <- SpatialFeaturePlot(
    object = control_obj,
    features = marker_panel_present,
    crop = FALSE,
    ncol = 2,
    alpha = c(0.1, 1)
  ) &
    manuscript_theme

  case_marker_spatial <- SpatialFeaturePlot(
    object = case_obj,
    features = marker_panel_present,
    crop = FALSE,
    ncol = 2,
    alpha = c(0.1, 1)
  ) &
    manuscript_theme

  save_plot(
    "plot_spatial_marker_expression_control.png",
    control_marker_spatial,
    width = FIG_TALL_WIDTH,
    height = FIG_TALL_HEIGHT
  )

  save_plot(
    "plot_spatial_marker_expression_case.png",
    case_marker_spatial,
    width = FIG_TALL_WIDTH,
    height = FIG_TALL_HEIGHT
  )
} else {
  message("None of the selected marker genes were found in the dataset.")
}


############################################################
# Final save
############################################################

# Save the final integrated Seurat object after clustering and marker plotting.
saveRDS(
  combined,
  file = file.path(out_dir, "final_integrated_spatial_object_no_DEG.rds")
)

message("Workflow complete.")
message("Outputs saved to: ", out_dir)
