#' Example multislice SpatialExperiment dataset
#'
#' A toy dataset derived from the DLPFC dataset provided by the
#' spatialLIBD package. The data have been subset and simplified
#' for demonstration purposes.
#'
#' The dataset consists of a list of SpatialExperiment objects,
#' each representing a tissue slice. Spot-level niche annotations
#' are stored in the \code{niche} column of \code{colData}.
#'
#' @format A list of SpatialExperiment objects
#' @source Derived from the DLPFC dataset in the spatialLIBD package
#' @references
#' Pardo B, Spangler A, Weber LM, Hicks SC, Jaffe AE, Martinowich K,
#' Maynard KR, Collado-Torres L (2022).
#' spatialLIBD: an R/Bioconductor package to visualize spatially-resolved
#' transcriptomics data. \emph{BMC Genomics}.
#' \doi{10.1186/s12864-022-08601-w}
#' @keywords datasets
"spe_list"
