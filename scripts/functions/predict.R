##################################################
# functions to  predict scRNAseq data
##################################################

# function to map project query on ref and make predictions based on Seurat integration
mapSeurat <- function(ref, query) {
    reference_list <- list(ref = ref, query = query)
    features <- Seurat::SelectIntegrationFeatures(object.list = reference_list)
    anchors <- Seurat::FindTransferAnchors(
        reference = reference_list$ref,
        query = reference_list$query,
        normalization.method = "LogNormalize",
        features = features
    )
    predictions <- Seurat::TransferData(
        anchorset = anchors,
        refdata = reference_list$ref$cluster
    )
    return(predictions)
}

# function to store predictions in seurat object
storePred <- function(predictions, label_col, score_col, seu_obj) {
    predictions_prep <-
        predictions |>
        tibble::rownames_to_column("barcode") |>
        dplyr::select(predicted.id, prediction.score.max, barcode) |>
        dplyr::mutate(
            predicted.id = ifelse(
                prediction.score.max < 0.3,
                "unknown",
                predicted.id
            )
        ) |>
        tibble::as_tibble() |>
        dplyr::rename(
            {{ label_col }} := predicted.id,
            {{ score_col }} := prediction.score.max
        )
    seu_obj@meta.data <-
        seu_obj@meta.data |>
        tibble::rownames_to_column("barcode") |>
        dplyr::left_join(predictions_prep, by = "barcode") |>
        tibble::column_to_rownames(var = "barcode")

    return(seu_obj)
}

# function to convert mouse rownames to human rownames
convertRownames <- function(seu_object) {
    lookup <- homologene::mouse2human(
        rownames(seu_object),
        db = homologene::homologeneData2
    )
    new_rownames <- lookup$humanGene[match(
        rownames(seu_object),
        lookup$mouseGene
    )]
    rownames(seu_object@assays$RNA@counts) <- new_rownames
    rownames(seu_object@assays$RNA@data) <- new_rownames
    # rownames(seu_object@assays$RNA@scale.data) <- new_rownames
    #remove columns with NA
    features_keep <- rownames(seu_object)[!is.na(rownames(seu_object))]
    obj_new <- subset(seu_object, features = features_keep)
    rownames(obj_new@assays$RNA@meta.features) <- rownames(obj_new)
    return(obj_new)
}

# function to predict scRNAseq data using Milbrandt et al. reference
predict_scRNAseq_milbrandt <- function(seu_obj) {
    seu_obj <- SeuratObject::JoinLayers(seu_obj)
    Seurat::DefaultAssay(seu_obj) <- "RNA"
    pns_sn_sciatic_milbrandt <- qs::qread(
        "/home/mischko/Documents/beruf/forschung/scRNA_reference/pns_atlas_milbrandt/pns_sn_sciatic_GSE182098.qs",
        nthreads = 4
    )
    human_pns_sciatic_milbrandt <- convertRownames(pns_sn_sciatic_milbrandt)
    predictions_milbrandt <- mapSeurat(
        ref = human_pns_sciatic_milbrandt,
        query = seu_obj
    )
    seu_obj <- storePred(
        predictions_milbrandt,
        label_col = "milbrandt_sciatic_label",
        score_col = "milbrandt_sciatic_score",
        seu_obj = seu_obj
    )
    seu_obj
}

predict_scRNAseq_heming <- function(seu_obj) {
    seu_obj <- SeuratObject::JoinLayers(seu_obj)
    Seurat::DefaultAssay(seu_obj) <- "RNA"
    pns_sn_heming <- qs::qread(
        "/home/mischko/Documents/beruf/forschung/seed/sn_sural_2023_10/objects/sc_merge_small.qs",
        nthreads = 4
    )
    predictions_heming <- mapSeurat(
        ref = pns_sn_heming,
        query = seu_obj
    )
    seu_obj <- storePred(
        predictions_heming,
        label_col = "heming_label",
        score_col = "heming_score",
        seu_obj = seu_obj
    )
    seu_obj
}
