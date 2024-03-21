#' @title MACE consensus
#'
#' @description This function performs the malignant cell (MACE) annotation, identifying each cell being normal, normal Epithelial or tumor based on the consensus method including both the Reference-based scATOMIC and the CNV-based SCEVAN methods.
#'
#' @param Obj A Seurat object containing the input single-cell data with it's corresponding metadata.
#' @param donor_id_col A string identifying the metadata column name where the sample identification is stored, this is important because the method is performed sample by sample.
#' @param Samples_to_annotate A list of unique sample ids to perform the cancer annotation for, has to be found in the donor_id_col. if NULL (default), the annotation will be performed for all samples.
#'
#' @return A Seurat Object containing the original input in addition to REF_cancer_annot and Consensus_cancer_annot metadata columns containing the MACE annotations.
#'
#' @examples
#' \dontrun{
#' Object <- MACE_consensus(Obj=Object, donor_id_col = "patient_id")
#' }
#'
#' @export
#'

MACE_consensus <- function(Obj,donor_id_col, Samples_to_annotate=NULL){
  # Run MACE_REF first, if it does not already exist
  if(is.element("REF_cancer_annot",colnames(Obj@meta.data))){
    print("Skipping REF_based annotation, as it already exists")
  }
  else{
    Obj <- MACE_REF(Obj, donor_id_col, Samples_to_annotate)
  }

  # set Idents of obj to donor/sample id column
  Idents(Obj) <- Obj@meta.data[donor_id_col]

  # create annotation column for the entire object
  df_CNV_annot = data.frame(matrix(nrow= dim(Obj)[2], ncol= 1))
  colnames(df_CNV_annot) <- 'CNVwnormal_cancer_annot'
  rownames(df_CNV_annot) <- colnames(Obj)

  # check if CNV_cancer_annot already exist
  if(is.element("CNVwnormal_cancer_annot",colnames(Obj@meta.data))){
    df_CNV_annot[,"CNVwnormal_cancer_annot"] = Obj@meta.data$CNVwnormal_cancer_annot
  }

  # Decide on samples list
  if(is.null(Samples_to_annotate)){
    Samples_to_annotate = unique(Idents(Obj))
  }

  # RUN SCEVAN annotation per sample using normal annotation from scATOMIC (REF)
  for (id in Samples_to_annotate){
    print(id)
    Obj.downsampled <- subset(Obj, idents= id)

    if (sum(Obj.downsampled$REF_cancer_annot == 'normal')>0){
      norm_cell_names = colnames(Obj.downsampled[,Obj.downsampled$REF_cancer_annot == 'normal'])
      if (length(norm_cell_names) > 50){
        norm_cell_names_to_use = norm_cell_names[sample(1:length(norm_cell_names),50,replace=F)]
        norm_cell_names_to_exclude = setdiff(norm_cell_names,norm_cell_names_to_use)
        Obj.downsampled = Obj.downsampled[,!(colnames(Obj.downsampled) %in% norm_cell_names_to_exclude)]
      }
      else{
        norm_cell_names_to_use = norm_cell_names
      }
    }
    else{
      norm_cell_names_to_use = NULL
    }

    results <- SCEVAN::pipelineCNA(Obj.downsampled@assays$RNA@counts, sample = id, par_cores = 20, norm_cell= norm_cell_names_to_use, SUBCLONES = FALSE, plotTree = FALSE)
    df_CNV_annot[rownames(results),"CNVwnormal_cancer_annot"] <- results$class
  }

  Obj <- AddMetaData(Obj, metadata=df_CNV_annot$CNVwnormal_cancer_annot, col.name="CNVwnormal_cancer_annot")

  # Generate consensus annotation
  df_concensus_annot = data.frame(matrix(nrow= dim(Obj)[2], ncol= 1))
  colnames(df_concensus_annot) <- 'Consensus_cancer_annot'
  rownames(df_concensus_annot) <- colnames(Obj)

  df_concensus_annot[Obj$REF_cancer_annot == 'filtered',"Consensus_cancer_annot"] <- 'filtered'
  df_concensus_annot[Obj$REF_cancer_annot == 'normal',"Consensus_cancer_annot"] <- 'normal'

  Obj$CNVwnormal_cancer_annot[is.na(Obj$CNVwnormal_cancer_annot == 'filtered' & !Obj$REF_cancer_annot == 'normal')] <- 'filtered'  # non normal cells that has no CNV annotation

  df_concensus_annot[(Obj$CNVwnormal_cancer_annot == 'filtered' & !Obj$REF_cancer_annot == 'normal'),"Consensus_cancer_annot"] <- 'filtered'
  df_concensus_annot[(Obj$CNVwnormal_cancer_annot == 'normal' & Obj$REF_cancer_annot == 'unclassified'),"Consensus_cancer_annot"] <- 'normal'
  df_concensus_annot[(Obj$CNVwnormal_cancer_annot == 'normal' & Obj$REF_cancer_annot == 'tumor'),"Consensus_cancer_annot"] <- 'normal Epithelial'
  df_concensus_annot[(Obj$CNVwnormal_cancer_annot == 'tumor' & !Obj$REF_cancer_annot == 'filtered' & !Obj$REF_cancer_annot == 'normal'),"Consensus_cancer_annot"] <- 'tumor'

  # Delete intermediate CNV annotation
  Obj$CNVwnormal_cancer_annot <- NULL

  # Add metadata column and return Seurat object
  Obj <- AddMetaData(Obj, metadata=df_concensus_annot$Consensus_cancer_annot, col.name="Consensus_cancer_annot")
  return(Obj)
}