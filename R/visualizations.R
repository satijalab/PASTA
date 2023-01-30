

#' Coverage Plot for genes based on residuals
#'
#' Calculate avergae length of tandem 3'UTRs 
#' @param object Seurat object
#' @param assay Name of polyAsite.Assay 
#' @param gene Name of gene to plot. Will extract all 
#' @param gene.names Column containing the gene where each polyA site is annotated.
#' @param extend.upstream 
#' @param extend.downstream
#' @importFrom utils head tail
#' @export
#' @concept visualization
#' 

PolyACoveragePlot <- function(object, 
                           assay = "polyA", 
                           gene = NULL,
                           gene.names = "symbol", 
                           extend.downstream = 500, 
                           extend.upstream = 500, 
                           annotation = "transcript", 
                           ...)
{
  meta <- object[[assay]][[]]
  meta.sub <- meta[rownames(object[[assay]]@scale.data),] #use scale data 
  tmp <- data.frame(peak = rownames(subset(meta.sub, get(gene.names) == gene)))
  tmp$start <- vapply(strsplit(tmp$peak, "-", fixed = TRUE), "[", "", 2)
  tmp$end <- vapply(strsplit(tmp$peak, "-", fixed = TRUE), "[", "", 3)
  tmp$chr <- vapply(strsplit(tmp$peak, "-", fixed = TRUE), "[", "", 1)
  tmp <- tmp[order(tmp$end),]
  range <- paste0(head(tmp$chr, 1), "-", head(tmp$star, 1), "-", tail(tmp$end, 1))
  CoveragePlot(object, region = range, assay = assay, extend.downstream = extend.downstream, 
               extend.upstream = extend.upstream, annotation = annotation, ...)
  
}

