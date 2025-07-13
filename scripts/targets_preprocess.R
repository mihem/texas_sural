# Set target options
tar_option_set(
  packages = c(
    "Seurat",
    "tidyverse",
    "scDblFinder",
    "scMisc",
    "tidyverse",
    "readxl",
    "writexl",
    "hdf5r"
  ),
  format = "qs"
)

# Define the targets pipeline
targets_preprocess <-
  tar_plan(
    sc_list = create_seurat_objects(),
    qc_plots = plot_qc(sc_list),
    sc_list_doublets = detect_doublets(sc_list),
    doublet_plots = plot_doublets(sc_list_doublets),
    sc_filter = filter_cells(sc_list_doublets)
    # filtered_plots = plot_filtered_cells(sc_filter, "results/qc"),
    # sc_merge_pre = merge_seurat_objects(sc_filter),
    # sc_merge = normalize_data(sc_merge_pre)
  )
