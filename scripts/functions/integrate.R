##################################################
# functions to integrate data
##################################################

# stancas integration semi-supervised based on azimuth annotation all
integrate_stacas_ss <- function(seu_obj) {
    stacas_obj <- Seurat::SplitObject(
        seu_obj,
        split.by = "sample"
    )
    stacas_obj <- STACAS::Run.STACAS(
        stacas_obj,
        cell.labels = "heming_label"
    )
    seu_obj$stacas <- stacas_obj$pca
    seu_obj
}

# umap ----
run_umap <- function(seu_obj) {
    seu_obj <- Seurat::RunUMAP(
        seu_obj,
        reduction = "stacas",
        dims = 1:30
    )
}

# add information to Seurat object
add_info <- function(seu_obj) {
    heming_small_misc <- qs::qread(
        "/home/mischko/Documents/beruf/forschung/seed/sn_sural_2023_10/objects/sc_merge_small_misc.qs"
    )
    seu_obj@misc$heming_cluster_col <- heming_small_misc$cluster_col
    seu_obj@misc$heming_cluster_order <- heming_small_misc$cluster_order
    seu_obj$level2 <- "DPN_Texas"
    seu_obj$cluster <- seu_obj$heming_label
    seu_obj
}

# plot umaps prediction
plot_umap <- function(seu_obj, group_by) {
    dir.create(file.path("results", "umap"), showWarnings = FALSE)
    umap <- Seurat::DimPlot(
        seu_obj,
        reduction = "umap",
        group.by = group_by,
        raster = FALSE,
        pt.size = .1,
        alpha = .1,
        cols = seu_obj@misc$heming_cluster_col,
        label = TRUE
    ) +
        Seurat::NoLegend() +
        ggplot2::ggtitle(group_by) +
        scMisc::theme_rect()
    ggplot2::ggsave(
        file.path(
            "results",
            "umap",
            paste0("umap_", group_by, ".png")
        ),
        plot = umap,
        width = 20,
        height = 8
    )
    return(umap)
}

# Combined texas and heming Seurat objects
merge_texas_heming <- function(seu_obj) {
    pns_sn_heming <- qs::qread(
        "/home/mischko/Documents/beruf/forschung/seed/sn_sural_2023_10/objects/sc_merge.qs",
        nthreads = 4
    )

    # merge heming and texas Seurat objects
    seu_obj_merge <- merge(
        x = pns_sn_heming,
        y = seu_obj,
        merge.data = FALSE
    )

    seu_obj_merge <-
        Seurat::DietSeurat(
            seu_obj_merge,
            layers = c("data"),
            assays = "RNA",
            dimreducs = NULL
        )

    # reduce the size of the Seurat object by creating empty data slots
    seu_obj_merge$RNA$data.1 <-
        Matrix::Matrix(
            0,
            nrow = nrow(seu_obj_merge),
            ncol = ncol(seu_obj_merge),
            sparse = TRUE
        )

    seu_obj_merge$RNA$data.2 <-
        Matrix::Matrix(
            0,
            nrow = nrow(seu_obj_merge),
            ncol = ncol(seu_obj_merge),
            sparse = TRUE
        )

    seu_obj_merge@misc$heming_cluster_col <- seu_obj@misc$heming_cluster_col
    seu_obj_merge@misc$heming_cluster_order <- seu_obj@misc$heming_cluster_order
    return(seu_obj_merge)
}