##################################################
# targets pipeline to calculcate abundance
##################################################

# Set target options
tar_option_set(
    packages = c(
        "Seurat",
        "tidyverse",
        "scMisc"
    ),
    format = "qs"
)

# Define the targets pipeline
targets_abundance <-
    tar_plan(
        # merge
        # abundance = plot_abundance(sc_merge_umap)
    )
