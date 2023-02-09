
 #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 # Functions
 #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
#' RunPCA only PolyA site assay
#' @param object A polyAsiteAssay
#' @param features Features to run PCA on
#' @param do.scale Scale residuals before performing PCA. Default is TRUE.
#' @param do.center Center residuals before performing PCA. Default is TRUE.
#' @param residuals.max Clip residuals above this value before performing PCA (default is 10)
#' @param residuals.min Clip residuals below this value before performing PCA (default is -10)
#' @param verbose Verbose output
#'
#' @rdname RunPCA 
#' @concept dimensional_reduction
#' @export
#' @method RunPCA polyAsiteAssay
#'
RunPCA.polyAsiteAssay <- function(
  object,
  features = NULL,
  do.scale = TRUE,
  do.center = TRUE,
  residuals.max = 10,
  residuals.min = -10,
  verbose=TRUE
) {

  residual.matrix <- object@scale.data
  if (is.null(features) &  !is.null(VariableFeatures(object)) ) {
    features <- VariableFeatures(object)
  }
  if (is.null(features) &  is.null(VariableFeatures(object)) ) {
    features <- rownames(residual.matrix)
  }

  residual.matrix <- residual.matrix[features,]
  if (verbose & do.scale & do.center) {
  message("Scaling and centering residual matrix")
  #add conditional messaging
  }
  #scaling by features
  residual.matrix <- t(scale(t(residual.matrix), center=do.center, scale= do.scale ))


  if (!is.null(residuals.max)) {
    if (verbose) {
    message(paste0("Clipping residuals above ", residuals.max))
    }
    residual.matrix[residual.matrix > residuals.max] <- residuals.max
  }

  if (!is.null(residuals.min)) {
    if (verbose) {
      message(paste0("Clipping residuals below ", residuals.min))
    }
    residual.matrix[residual.matrix < residuals.min] <- residuals.min
  }

  #make a new temporary assay to run PCA
  pcs <- RunPCA(residual.matrix)
  return(pcs)
}