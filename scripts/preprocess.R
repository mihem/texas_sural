##################################################
# preprocess scRNAseq data
##################################################

# load libraries ------
library(tidyverse)
library(readxl)
library(Seurat)
library(presto)
library(writexl)
library(scMisc)
library(hdf5r)
library(qs2)
library(scDblFinder)

# optional libraries for better coding
library(languageserver)
library(httpgd)

# general settings  ----
future::plan("multicore", workers = 6)
options(future.globals.maxSize = 16000 * 1024^2)

# find  raw files and match to names from lookup table ---
h5_path <- list.files(
    file.path("raw"),
    pattern = ".h5",
    recursive = TRUE,
    full.names = TRUE
)

seq_names <-
    tibble(name = h5_path) |>
    dplyr::mutate(
        name = gsub(
            x = name,
            pattern = "raw/([^/]+)/.*",
            replacement = "\\1"
        )
    ) |>
    pull(name)

# read in data and create Seurat object -----
sc_list_matrix <-
    lapply(h5_path, scMisc::ReadCellBender_h5) |>
    setNames(seq_names)

sc_list <- lapply(sc_list_matrix, FUN = function(x) {
    CreateSeuratObject(
        counts = x,
        min.cells = 3,
        min.features = 200,
        project = "PNP"
    )
})

# plot QC: MT and genes ----
for (i in seq_along(sc_list)) {
    sc_list[[i]][["percent_mt"]] <- PercentageFeatureSet(
        sc_list[[i]],
        pattern = "^MT"
    )
}

plot1 <- vector("list", length = length(sc_list))

for (i in seq_along(sc_list)) {
    plot1[[i]] <- FeatureScatter(
        object = sc_list[[i]],
        feature1 = "nCount_RNA",
        feature2 = "percent_mt",
        raster = FALSE
    ) +
        labs(title = names(sc_list)[[i]]) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
        NoLegend()
}

plot1_patch <- patchwork::wrap_plots(plot1, ncol = 4)

ggsave(
    file.path("results", "qc", "mt.png"),
    width = 15,
    height = 5,
)

plot2 <- vector("list", length = length(sc_list))

for (i in seq_along(sc_list)) {
    plot2[[i]] <- FeatureScatter(
        object = sc_list[[i]],
        feature1 = "nCount_RNA",
        feature2 = "nFeature_RNA"
    ) +
        labs(title = names(sc_list)[[i]]) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
        NoLegend()
}

plot2_patch <- patchwork::wrap_plots(plot2, ncol = 4)

ggsave(
    file.path("results", "qc", "genes.png"),
    width = 15,
    height = 5,
)

#calculcate doubles and add scDblFinder scores to seurat objects
doubletFun <- function(seu_obj) {
    sce <- Seurat::as.SingleCellExperiment(seu_obj)
    sce <- scDblFinder::scDblFinder(sce)
    seu_obj$scDblFinder.score <- sce$scDblFinder.score
    seu_obj$scDblFinder.class <- sce$scDblFinder.class
    return(seu_obj)
}

sc_list <- lapply(sc_list, doubletFun)

#plot scDblFinder
plot3 <- vector("list")

for (i in seq_along(sc_list)) {
    plot3[[i]] <- FeatureScatter(
        object = sc_list[[i]],
        feature1 = "nCount_RNA",
        feature2 = "nFeature_RNA",
        group.by = "scDblFinder.class"
    ) +
        labs(title = names(sc_list)[i]) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
}

plot3_patch <- patchwork::wrap_plots(plot3, ncol = 4)

ggsave(
    file.path("results", "qc", "doublet.png"),
    width = 15,
    height = 5
)

#create table with doublet rate
doublet_tbl <- purrr::map_dfr(
    purrr::map(sc_list, "meta.data"),
    dplyr::count,
    scDblFinder.class
) |>
    dplyr::mutate(sample = rep(names(sc_list), each = 2)) |>
    tidyr::pivot_wider(names_from = "scDblFinder.class", values_from = "n") |>
    dplyr::mutate(doublet_pct = doublet / (singlet + doublet) * 100)

write_csv(doublet_tbl, file.path("results", "qc", "doublets.csv"))

# filter low quality cells and doublets ---
filter_df <- readr::read_csv(file.path("lookup", "filter_df.csv"))

sc_filter <- vector("list", length(sc_list))

for (i in seq_along(sc_list)) {
    sc_filter[[i]] <-
        subset(
            sc_list[[i]],
            subset = nFeature_RNA > 200 &
                nFeature_RNA < filter_df$rna[[i]] &
                percent_mt < filter_df$mt[[i]] &
                scDblFinder.class == "singlet"
        )
}

names(sc_filter) <- names(sc_list)

# check filter settings ---
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
    file.path("results", "qc", "mt_post.png"),
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
    file.path("results", "qc", "genes_post.png"),
    width = 15,
    height = 5
)

# merge seurat objects  ------
sc_merge_pre <- merge(
    x = sc_filter[[1]],
    y = sc_filter[-1],
    merge.data = TRUE,
    add.cell.ids = names(sc_list)
)

sc_merge_pre$sample <- gsub(x = colnames(sc_merge_pre), pattern = "(SN_Dallas_\\d+)_.*", replacement = "\\1")

# this needs to be rejoined and then split again to get the correct names (better naming and required for integration)
sc_merge_pre <- JoinLayers(sc_merge_pre)
sc_merge_pre <- split(x = sc_merge_pre, f = sc_merge_pre$sample)

# # add metadata
# seq_metadata <-
#     lookup |>
#     dplyr::select(
#         patient,
#         sex,
#         age,
#         group,
#         diagnosis,
#         incat_at_lumbar_puncture,
#         incat_follow_up,
#         onls_at_lumbar_puncture,
#         onls_follow_up,
#         mrc_sum_score_60_at_lumbar_puncture,
#         mrc_sum_score_60_follow_up,
#         icu
#     )

# seq_metadata <-
#     seq_metadata |>
#     dplyr::select(patient, batch)

# sc_merge_pre@meta.data <-
#     sc_merge_pre@meta.data |>
#     tibble::rownames_to_column("barcode") |>
#     dplyr::left_join(seq_metadata, by = "patient") |>
#     tibble::column_to_rownames(var = "barcode")

# sc_merge_pre$batch <- paste0(sc_merge_pre$batch, "_", sc_merge_pre$tissue)

# str(sc_merge@meta.data)
# str(sc_merge_pre@meta.data)

# qc metrics  -----
metrics_files <- list.files(
    file.path("raw"),
    pattern = "metrics_summary.csv",
    recursive = TRUE,
    full.names = TRUE
)

metrics_data <-
    purrr::map_df(metrics_files, read_csv) |>
    mutate(sample = seq_names, .before = `Estimated Number of Cells`) |>
    arrange(sample)

write_csv(metrics_data, file.path("results", "qc", "cellranger_metrics.csv"))

count_cells <-
    purrr::map_df(sc_list, dim) |>
    dplyr::slice(2) |>
    tidyr::pivot_longer(everything(), names_to = "sample") |>
    dplyr::left_join(dplyr::count(sc_merge_pre@meta.data, sample)) |>
    dplyr::rename(before = value, after = n)

write_csv(count_cells, file.path("results", "qc", "count_cells.csv"))

count_genes <-
    data.frame(
        feature = sc_merge_pre@meta.data$nFeature_RNA,
        sample = sc_merge_pre@meta.data$sample
    ) |>
    dplyr::group_by(sample) |>
    dplyr::summarize(median_genes_after = median(feature))

write_csv(count_genes, file.path("results", "qc", "count_genes.csv"))

# normalize ----
sc_merge <- NormalizeData(
    sc_merge_pre,
    verbose = TRUE,
    normalization.method = "LogNormalize",
    scale.factor = 10000
)

sc_merge <-
    sc_merge |>
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) |>
    ScaleData() |>
    RunPCA()

