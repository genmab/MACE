#' @title MACE REF
#'
#' @description This function performs the malignant cell (MACE) annotation, identifying each cell being normal or tumor based on the Reference-based method scATOMIC.
#'
#' @param Obj A Seurat object containing the input single-cell data with it's corresponding metadata.
#' @param donor_id_col A string identifying the metadata column name where the sample identification is stored, this is important because the method is performed sample by sample.
#' @param Samples_to_annotate A list of unique sample ids to perform the cancer annotation for, has to be found in the donor_id_col. if NULL (default), the annotation will be performed for all samples.
#'
#' @return A Seurat Object containing the original input in addition to REF_cancer_annot metadata column containing the MACE annotations.
#'
#' @examples
#' \dontrun{
#' Object <- MACE_REF(Obj=Object, donor_id_col = "patient_id")
#' }
#'
#' @export
#' @import Seurat scATOMIC
#'

MACE_REF <- function(Obj,donor_id_col, Samples_to_annotate=NULL){
  # set Idents of obj to donor/sample id column
  Idents(Obj) <- Obj@meta.data[donor_id_col]

  # create annotation column for the entire object
  df_REF_annot = data.frame(matrix(nrow= dim(Obj)[2], ncol= 1))
  colnames(df_REF_annot) <- 'REF_cancer_annot'
  rownames(df_REF_annot) <- colnames(Obj)

  # check if REF_cancer_annot already exist
  if(is.element("REF_cancer_annot",colnames(Obj@meta.data))){
    df_REF_annot[,"REF_cancer_annot"] = Obj@meta.data$REF_cancer_annot
  }

  # Decide on samples list
  if(is.null(Samples_to_annotate)){
    Samples_to_annotate = unique(Idents(Obj))
  }

  # RUN SCEVAN annotation per sample
  for (id in Samples_to_annotate){
    print(id)
    Obj.downsampled <- subset(Obj, idents= id)

    count_data <- Obj.downsampled@assays$RNA@counts
    pct_mt <- colSums(as.matrix(count_data[grep("^MT-", row.names(count_data)),]))/colSums(as.matrix(count_data)) * 100
    nFeatureRNA <- colSums(as.matrix(count_data > 0))
    count_data <- count_data[, names(which(pct_mt < 25))]
    count_data <- count_data[, intersect(names(which(nFeatureRNA > 200)), colnames(count_data))]

    cell_predictions <- scATOMIC::run_scATOMIC(count_data, fine_grained_T = F)
    results <- create_summary_matrix(prediction_list = cell_predictions, use_CNVs = F, modify_results = F, raw_counts = count_data, min_prop = 0.5, fine_grained_T = F)

    Pred_class = data.frame(matrix(ncol=1,nrow=dim(Obj.downsampled)[2]))
    colnames(Pred_class) <- c("Pred.class")
    rownames(Pred_class) <- colnames(Obj.downsampled)

    Pred_class$Pred.class <- 'filtered'
    Pred_class[rownames(results),'Pred.class'] <- results$layer_2
    Pred_class$Pred.class[Pred_class$Pred.class == 'unclassified_any_cell' | Pred_class$Pred.class == 'unclassified_normal_or_cancer_tissue'] <- 'unclassified'
    Pred_class$Pred.class[Pred_class$Pred.class == 'Non Stromal Cell'] <- 'tumor'
    Pred_class$Pred.class[!(Pred_class$Pred.class %in% c('filtered','unclassified','tumor'))] <- 'normal'

    df_REF_annot[rownames(Pred_class),"REF_cancer_annot"] <- Pred_class$Pred.class
  }

  # Add metadata column and return Seurat object
  Obj <- AddMetaData(Obj, metadata=df_REF_annot$REF_cancer_annot, col.name="REF_cancer_annot")
  return(Obj)
}