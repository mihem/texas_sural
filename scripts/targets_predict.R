##################################################
# targets pipeline to predict scRNAseq data
##################################################

# Set target options
tar_option_set(
    packages = c(
        "Seurat",
        "BPCells",
        "tidyverse",
        "qs"
    ),
    format = "qs"
)

# Define the targets pipeline
targets_predict <-
    tar_plan(
        sc_merge_predict_milbrandt = predict_scRNAseq_milbrandt(sc_merge),
        sc_merge_predict_heming = predict_scRNAseq_heming(sc_merge_predict_milbrandt)
    )
