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

# # umap ----
run_umap <- function(seu_obj) {
    seu_obj <- Seurat::RunUMAP(
        seu_obj,
        reduction = "stacas",
        dims = 1:30
    )
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
        # cols = colors_dutch,
        label = TRUE
    ) +
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
}
