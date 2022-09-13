
#' Find Markers for a polyA site Assay
#' 
#' @rdname FindMarkers
#' @concept differential_expression
#' @export
#' @method FindMarkers polyAsiteAssay


# FindMarkers.polyAsiteAssay <- function(
#     object, 
#     slot = "scale.data", 
#     cells.1 = NULL,
#     cells.2 = NULL,
#     features = NULL, 
#     ...) {
#   object <- NextMethod()
#   m <- FindMarkers(object = object, slot=slot, cells.1 = cells.1, cells.2 = cells.2, features = features, ...)
#   m$symbol <- object[["symbol"]][rownames(m),]
#   return(m)
# }


