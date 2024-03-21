#' @title MACE CNV
#'
#' @description This function performs the malignant cell (MACE) annotation, identifying each cell being normal or tumor based on the CNV-based method SCEVAN.
#'
#' @param Obj A Seurat object containing the input single-cell data with it's corresponding metadata.
#' @param donor_id_col A string identifying the metadata column name where the sample identification is stored, this is important because the method is performed sample by sample.
#' @param Samples_to_annotate A list of unique sample ids to perform the cancer annotation for, has to be found in the donor_id_col. if NULL (default), the annotation will be performed for all samples.
#'
#' @return A Seurat Object containing the original input in addition to CNV_cancer_annot metadata column containing the MACE annotations.
#'
#' @examples
#' \dontrun{
#' Object <- MACE_CNV(Obj=Object, donor_id_col = "patient_id")
#' }
#'
#' @export
#' @import SCEVAN
#'

MACE_CNV <- function(Obj,donor_id_col, Samples_to_annotate=NULL){
  # set Idents of obj to donor/sample id column
  Idents(Obj) <- Obj@meta.data[donor_id_col]

  # create annotation column for the entire object
  df_CNV_annot = data.frame(matrix(nrow= dim(Obj)[2], ncol= 1))
  colnames(df_CNV_annot) <- 'CNV_cancer_annot'
  rownames(df_CNV_annot) <- colnames(Obj)

  # check if CNV_cancer_annot already exist
  if(is.element("CNV_cancer_annot",colnames(Obj@meta.data))){
    df_CNV_annot[,"CNV_cancer_annot"] = Obj@meta.data$CNV_cancer_annot
  }

  # Decide on samples list
  if(is.null(Samples_to_annotate)){
    Samples_to_annotate = unique(Idents(Obj))
  }

  # RUN SCEVAN annotation per sample
  for (id in Samples_to_annotate){
    print(id)
    Obj.downsampled <- subset(Obj, idents= id)

    results <- SCEVAN::pipelineCNA(Obj.downsampled@assays$RNA@counts, sample = id, par_cores = 20, SUBCLONES = FALSE, plotTree = FALSE)
    df_CNV_annot[rownames(results),"CNV_cancer_annot"] <- results$class
  }

  # Add metadata column and return Seurat object
  Obj <- AddMetaData(Obj, metadata=df_CNV_annot$CNV_cancer_annot, col.name="CNV_cancer_annot")
  return(Obj)
}