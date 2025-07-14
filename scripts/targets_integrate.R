##################################################
# targets pipeline to integrate scRNAseq data
##################################################

# Set target options
tar_option_set(
    packages = c(
        "Seurat",
        "STACAS",
        "tidyverse",
        "scMisc"
    ),
    format = "qs"
)

# Define the targets pipeline
targets_integrate <-
    tar_plan(
        sc_merge_stacas = integrate_stacas_ss(sc_merge_predict_heming),
        sc_merge_heming_colors = add_heming_small_misc_colors(sc_merge_stacas),
        sc_merge_umap = run_umap(sc_merge_heming_colors),
        umap_heming = plot_umap(sc_merge_umap, "heming_label")
    )
