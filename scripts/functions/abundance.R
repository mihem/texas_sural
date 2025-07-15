##################################################
# functions to calculate abundance
##################################################
plot_abundance <- function(seu_obj) {
    dir.create(file.path("results", "abundance"), showWarnings = FALSE)
    abundance_plot <- scMisc::stackedPlot(
        object = seu_obj,
        x_axis = "sample",
        y_axis = "heming_label",
        x_order = unique(seu_obj$sample),
        y_order = seu_obj@misc$heming_cluster_order,
        color = seu_obj@misc$heming_cluster_col,
        width = 10
    )
    ggplot2::ggsave(
        file.path("results", "abundance", "abundance_plot.png"),
        plot = abundance_plot,
        width = 5,
        height = 8
    )
    return(abundance_plot)
}

# abundance_plot <- scMisc::stackedPlot(
#     object = seu_obj,
#     x_axis = "sample",
#     y_axis = "cluster",
#     x_order = unique(sc_merge_heming_texas_small$sample),
#     y_order = sc_merge_heming_texas_small@misc$heming_cluster_order,
#     color = sc_merge_heming_texas_small@misc$heming_cluster_col,
#     width = 10
# )

# sample_lookup <-
#     readr::read_csv(
#         "/home/mischko/Documents/beruf/forschung/seed/sn_sural_2023_10/lookup/sample_lookup.csv"
#     ) |>
#     mutate(level0 = if_else(level1 == "CTRL", "CTRL", "PNP"))
# propeller_PNP_CTRL <-
#     scMisc::propellerCalc(
#         seu_obj1 = sc_merge_heming_texas,
#         condition1 = "DPN_Texas",
#         condition2 = "CTRL",
#         cluster_col = "cluster",
#         meta_col = "level2",
#         lookup = sample_lookup,
#         sample_col = "sample",
#         formula = "~0 + level2",
#         min_cells = 30
#     )

# scMisc::plotPropeller(
#     data = propeller_PNP_CTRL,
#     color = sc_merge@misc$cluster_col,
#     filename = "PNP_CTRL",
#     FDR = 0.1
# )

# scMisc::dotplotPropeller(
#     data = propeller_PNP_CTRL,
#     color = sc_merge@misc$cluster_col,
#     filename = "PNP_CTRL",
# )
