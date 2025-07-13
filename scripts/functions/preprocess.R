##################################################
# preprocess scRNAseq data
##################################################

# Function to create Seurat objects and calculate QC metrics
create_seurat_objects <- function() {
  h5_path <- list.files(
    file.path("raw"),
    pattern = ".h5",
    recursive = TRUE,
    full.names = TRUE
  )
  seq_names <- data.frame(name = h5_path) |>
    dplyr::mutate(
      name = gsub(
        x = name,
        pattern = "raw/([^/]+)/.*",
        replacement = "\\1"
      )
    ) |>
    dplyr::pull(name)
  sc_list_matrix <- lapply(h5_path, scMisc::ReadCellBender_h5) |>
    setNames(seq_names)
  sc_list <- lapply(sc_list_matrix, FUN = function(x) {
    CreateSeuratObject(
      counts = x,
      min.cells = 3,
      min.features = 200,
      project = "PNP"
    )
  })
  for (i in seq_along(sc_list)) {
    sc_list[[i]][["percent_mt"]] <- PercentageFeatureSet(
      sc_list[[i]],
      pattern = "^MT"
    )
  }
  sc_list
}

# Function to detect doublets
detect_doublets <- function(seu_obj_list) {
  doubletFun <- function(seu_obj) {
    sce <- Seurat::as.SingleCellExperiment(seu_obj)
    sce <- scDblFinder::scDblFinder(sce)
    seu_obj$scDblFinder.score <- sce$scDblFinder.score
    seu_obj$scDblFinder.class <- sce$scDblFinder.class
    seu_obj
  }
  lapply(seu_obj_list, doubletFun)
}

# Function to filter cells
filter_cells <- function(sc_list, filter_df) {
  sc_filter <- vector("list", length(sc_list))
  for (i in seq_along(sc_list)) {
    sc_filter[[i]] <- subset(
      sc_list[[i]],
      subset = nFeature_RNA > 200 &
        nFeature_RNA < filter_df$rna[[i]] &
        percent_mt < filter_df$mt[[i]] &
        scDblFinder.class == "singlet"
    )
  }
  names(sc_filter) <- names(sc_list)
  sc_filter
}

# Function to merge Seurat objects
merge_seurat_objects <- function(sc_filter) {
  sc_merge_pre <- merge(
    x = sc_filter[[1]],
    y = sc_filter[-1],
    merge.data = TRUE,
    add.cell.ids = names(sc_filter)
  )
  sc_merge_pre$sample <- gsub(
    x = colnames(sc_merge_pre),
    pattern = "(SN_Dallas_[0-9]+)_.*",
    replacement = "\1"
  )
  sc_merge_pre
}

# Function to normalize data
normalize_data <- function(sc_merge_pre) {
  sc_merge <- NormalizeData(
    sc_merge_pre,
    verbose = TRUE,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  sc_merge |>
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) |>
    ScaleData() |>
    RunPCA()
}

# Function to plot QC metrics
plot_qc <- function(seu_obj_list) {

  plot1 <- vector("list", length = length(seu_obj_list))
  for (i in seq_along(seu_obj_list)) {
    plot1[[i]] <- FeatureScatter(
      object = seu_obj_list[[i]],
      feature1 = "nCount_RNA",
      feature2 = "percent_mt",
      raster = FALSE
    ) +
      labs(title = names(seu_obj_list)[[i]]) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
      NoLegend()
  }
  plot1_patch <- patchwork::wrap_plots(plot1, ncol = 4)
  ggsave(
    file.path("results", "qc", "mt.png"),
    plot = plot1_patch,
    width = 15,
    height = 5
  )

  plot2 <- vector("list", length = length(seu_obj_list))
  for (i in seq_along(seu_obj_list)) {
    plot2[[i]] <- FeatureScatter(
      object = seu_obj_list[[i]],
      feature1 = "nCount_RNA",
      feature2 = "nFeature_RNA"
    ) +
      labs(title = names(seu_obj_list)[[i]]) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
      NoLegend()
  }
  plot2_patch <- patchwork::wrap_plots(plot2, ncol = 4)
  ggsave(
    file.path("results", "qc", "genes.png"),
    plot = plot2_patch,
    width = 15,
    height = 5
  )
}

# Function to plot doublets
plot_doublets <- function(seu_obj_list) {
  plot3 <- vector("list", length = length(seu_obj_list))
  for (i in seq_along(seu_obj_list)) {
    plot3[[i]] <- FeatureScatter(
      object = seu_obj_list[[i]],
      feature1 = "nCount_RNA",
      feature2 = "nFeature_RNA",
      group.by = "scDblFinder.class"
    ) +
      labs(title = names(seu_obj_list)[[i]]) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
  }
  plot3_patch <- patchwork::wrap_plots(plot3, ncol = 4)
  ggsave(
    file.path("results", "qc", "doublet.png"),
    plot = plot3_patch,
    width = 15,
    height = 5
  )
}

# Function to plot filtered cells
plot_filtered_cells <- function(sc_filter, output_dir) {
  plot4 <- vector("list", length = length(sc_filter))
  for (i in seq_along(sc_filter)) {
    plot4[[i]] <- FeatureScatter(
      object = sc_filter[[i]],
      feature1 = "nCount_RNA",
      feature2 = "percent_mt",
      raster = FALSE
    ) +
      labs(title = names(sc_filter)[[i]]) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
      NoLegend()
  }
  plot4_patch <- patchwork::wrap_plots(plot4, ncol = 4)
  ggsave(
    file.path(output_dir, "mt_post.png"),
    plot = plot4_patch,
    width = 15,
    height = 5
  )

  plot5 <- vector("list", length = length(sc_filter))
  for (i in seq_along(sc_filter)) {
    plot5[[i]] <- FeatureScatter(
      object = sc_filter[[i]],
      feature1 = "nCount_RNA",
      feature2 = "nFeature_RNA"
    ) +
      labs(title = names(sc_filter)[[i]]) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
      NoLegend()
  }
  plot5_patch <- patchwork::wrap_plots(plot5, ncol = 4)
  ggsave(
    file.path(output_dir, "genes_post.png"),
    plot = plot5_patch,
    width = 15,
    height = 5
  )
}
