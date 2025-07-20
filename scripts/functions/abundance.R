##################################################
# functions to calculate abundance
##################################################
# function to plot abundance
plot_abundance <- function(sc_heming_texas) {
    dir.create(file.path("results", "abundance"), showWarnings = FALSE)
    sc_heming_texas <- subset(sc_heming_texas, cluster != "unknown")
    scMisc::stackedPlot(
        object = sc_heming_texas,
        x_axis = "sample",
        y_axis = "cluster",
        x_order = unique(sc_heming_texas$sample),
        y_order = sc_heming_texas@misc$heming_cluster_order,
        color = sc_heming_texas@misc$heming_cluster_col,
        width = 10,
        dir_output = file.path("results", "abundance"),
    )
}

# function to create a sample lookup table
create_sample_lookup <- function(texas_file, heming_file) {
    sample_lookup_texas <- readr::read_csv(texas_file)
    sample_lookup_heming <- readr::read_csv(heming_file)

    sample_lookup <- dplyr::bind_rows(
        sample_lookup_heming,
        sample_lookup_texas
    )
    return(sample_lookup)
}

# function to calculate propeller DPN_Texas vs CTRL
plot_propeller_DPN_CTRL <- function(
    seu_obj,
    sample_lookup,
    condition1,
    condition2
) {
    seu_obj <- subset(seu_obj, cluster != "unknown")
    propeller_data <-
        scMisc::propellerCalc(
            seu_obj1 = seu_obj,
            condition1 = condition1,
            condition2 = condition2,
            cluster_col = "cluster",
            meta_col = "level2",
            lookup = sample_lookup,
            sample_col = "sample",
            formula = "~0 + level2 + age + sex",
            min_cells = 10
        )

    scMisc::plotPropeller(
        data = propeller_data,
        color = seu_obj@misc$heming_cluster_col,
        dir_output = file.path("results", "abundance"),
        filename = paste0(
            condition1,
            "_vs_",
            condition2
        ),
        FDR = 0.1
    )
}
