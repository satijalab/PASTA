

#' Coverage Plot for genes based on polyA sites.
#'
#' @param object Seurat object
#' @param region Region to plot. If NULL (default), will show all polyA sites within the 
#' gene for which polyA residuals have been calculated. If region is not specified and 
#' polyA residuals have not been calculated, all polyA sites within a gene will be shown.
#' @param assay Name of polyAsite.Assay Default is "polyA"
#' @param gene Name of gene to plot.
#' @param gene.names Column containing the gene where each polyA site is annotated.
#' @param extend.downstream How far to extend coverage downstream. Default is 500 bp. 
#' @param extend.downstream How far to extend coverage upstream. Default is 500 bp. 
#' @param annotation Default is "transcript". Use "gene" if you only want to see one isoform. 
#' @importFrom utils head tail
#' @export
#' @concept visualization
#' @seealso \code{\link[Signac]{CoveragePlot}}  for additional parameters
#' 

PolyACoveragePlot <- function(object, 
                              region.plot = NULL,
                              assay = "polyA", 
                              gene = NULL,
                              gene.names = "Gene_Symbol", 
                              extend.downstream = 500, 
                              extend.upstream = 500, 
                              annotation = "transcript", 
                           ...)
{
  
  if (!(is.null(region.plot))) {
      region.plot.tmp <- region.plot
  } else { # no region specified
    meta <- object[[assay]][[]]
     if ( dim(object[[assay]]@scale.data)[[1]] == 0) { #no polyA residuals
       message("PolyA Residuals not calculated and no region specified, will show all polyA sites in gene.")
       meta.sub <- meta
    } else { #use polyA residuals
      meta.sub <- meta[rownames(object[[assay]]@scale.data),] #use scale data 
    }
    tmp <- data.frame(peak = rownames(subset(meta.sub, get(gene.names) == gene)))
    tmp$start <- vapply(strsplit(tmp$peak, "-", fixed = TRUE), "[", "", 2)
    tmp$end <- vapply(strsplit(tmp$peak, "-", fixed = TRUE), "[", "", 3)
    tmp$chr <- vapply(strsplit(tmp$peak, "-", fixed = TRUE), "[", "", 1)
    tmp <- tmp[order(tmp$end),]
    region.plot.tmp <- paste0(head(tmp$chr, 1), "-", head(tmp$star, 1), "-", tail(tmp$end, 1))
  } 
  
  CoveragePlot(object, region = region.plot.tmp, assay = assay, extend.downstream = extend.downstream, 
               extend.upstream = extend.upstream, annotation = annotation, ...)
  
}

