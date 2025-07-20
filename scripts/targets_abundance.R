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
        # Track CSV files as dependencies
        abundance_stacked_sample = plot_abundance(sc_heming_texas),
        sample_lookup_texas_file = file.path("lookup", "sample_lookup_texas.csv"),
        sample_lookup_heming_file = file.path("lookup", "sample_lookup_heming.csv"),
        sample_lookup = create_sample_lookup(
            texas_file = sample_lookup_texas_file,
            heming_file = sample_lookup_heming_file
        ),
        propeller_DPN_CTRL = plot_propeller_DPN_CTRL(
            seu_obj = sc_heming_texas,
            sample_lookup = sample_lookup,
            condition1 = "DPN_Texas",
            condition2 = "CTRL"
        ),
        propeller_DPN_CIDP = plot_propeller_DPN_CTRL(
            seu_obj = sc_heming_texas,
            sample_lookup = sample_lookup,
            condition1 = "DPN_Texas",
            condition2 = "CIDP"
        ),
        propeller_DPN_CIAP = plot_propeller_DPN_CTRL(
            seu_obj = sc_heming_texas,
            sample_lookup = sample_lookup,
            condition1 = "DPN_Texas",
            condition2 = "CIAP"
        )
    )
