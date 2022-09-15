
#' Features to GRanges
#'
#' Transform List of Features to GRanges object
#' @param regions Vector of genomic coordinate ranges in format chr:start,end:strand
#' @param sep Vector of separators to use for genomic ranges string. First element separates chromosome
#' and start position, second element separates start position and end end position,
#' third element separates end position and strand.
#'
#' @return Returns a genomic ranges
#'
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame strand
#' @importFrom tidyr separate
#' @export
#' @concept utilities
FeaturesToGRanges <- function(regions, sep = c(":",",",":")) {
  ranges.df <- data.frame(ranges = regions)

  ranges.df <- separate(
    data = ranges.df,
    col = "ranges",
    sep = paste0(sep, collapse = "|"),
    into = c("chr", "start", "end","strand")
  )
  granges <- makeGRangesFromDataFrame(df = ranges.df)
  granges$str <- strand(granges)
  return(granges)
}

#' GRanges to String
#'
#' Converts GRanges object to vector of strings.
#' @param grange A GRanges object containing genomic coordinates and strand
#' @param sep   Vector of separators to use for genomic ranges string. First element separates chromosome
#' and start position, second element separates start position and end end position,
#' third element separates end position and strand
#'
#' @return Returns a character vector
#' 
#' @importFrom GenomicRanges makeGRangesFromDataFrame strand start end
#' 
#' @export
#' @concept utilities
GRangesToString <- function(grange, sep = c(":", ",", ":")) {
  regions <- paste0(
    as.character(x = seqnames(x = grange)),
    sep[[1]],
    start(x = grange),
    sep[[2]],
    end(x = grange),
    sep[[3]],
    strand(x = grange)
  )
  return(regions)
}


#' Calculate Tandem UTR Length Wrapper
#'
#' Calculate avergae length of tandem 3'UTRs 
#' @param object Seurat object
#' @param assay Name of polyAsite.Assay 
#' @param genes Vector of Gene names containing tandem UTRs. If Null, all annotated tandem UTRs are used (default). 
#' @param pseudocount Pseudocount to add for tandem UTR length calculation (default 0). 
#' @param min.reads.per.cell Use genes with minimum average read count (default 1). 
#'
#' @return Returns a matrix with average tandem UTR length
#' @import pbapply
#' @export
#' @concept utilities

CalculateTandemUTRlength <- function(object = NULL, assay = NULL, genes = NULL, pseudocount = 0, min.reads.per.cell = 1){
  if( !assay %in% Assays(object) ){
    stop(paste0(assay," assay is not present in object"))
  }
  if( !class(object[[assay]]) == "polyAsiteAssay"){
    stop(paste0(assay," assay is not a polyAsiteAssay"))
  }
  if( !assay == DefaultAssay(object)){
    DefaultAssay(object) <- assay
    message(paste0("Setting default assay to ",assay))
  }
  if( is.null(genes) ){
    # define tandem UTR genes
    idx = which( object[[assay]]@meta.features$is.tandem == TRUE )
    genes = object[[assay]]@meta.features$symbol[idx] %>% unique()
    # ! some sites are annotated to >1 gene.
    # Needs to be resolved later.
    # Ignored for now.
  }
  if( length(grep('is.tandem',colnames(object[[assay]]@meta.features))) == 0 ){
    stop("Exiting. Run AnnotatePASfromGTF() first.")
  }
  message(paste0("Looking for ",length(genes)," tandem UTR genes"))
  # get total raw counts per gene over cells 
  total <- rowSums( GetAssayData(object = object, assay = assay, slot = "data") )
  df <- data.frame(object[[assay]]@meta.features$symbol , 
                   object[[assay]]@meta.features$is.tandem ,
                   object[[assay]]@meta.features$tandem_utr_reference_length ,
                   object[[assay]]@meta.features$TandemUTRfraction ,
                   total)
  colnames(df) = c("genes","is.tandem","tandem_utr_reference_length","TandemUTRfraction","total")
  # only consider defined gene set
  df = df[which(df$genes %in% genes),]
  # only consider tandem UTR sites
  df = subset(df , is.tandem == TRUE)
  # remove sites with double gene annotation or single cleavage site in tandem UTR
  # ^^ see above, annotation needs to be resolved
  del=which( sapply(df$genes, FUN = function(j){strsplit(j,",")[[1]] %>% length()}) != 1)
  if(length(del)>0){
    message(paste0("Removed ",length(del)," sites with ambigous annotation"))
    df = df[-del,]
  }
  del=which( table(df$genes) == 1 )
  if(length(del)>0){
    message(paste0("Removed ",length(del)," genes with <2 sites per tandem UTR after filtering"))
    df = subset( df , ! genes %in% names(del))
  }
  
  L = pbsapply( unique(df$genes) ,function(gene) calcUTRlength(gene,df,object = object, assay = assay, pc = pseudocount , Min.Reads.Per.Cell = min.reads.per.cell))
  UTRlength.matrix = do.call( cbind , L[which( sapply(L,is.null) == FALSE )] )
  
  message(paste0("Calculated tandem UTR length for ", ncol(UTRlength.matrix),"/",length(unique(df$genes))," filtered tandem UTR genes"))
  
  return(UTRlength.matrix)
}

#' Calculate Tandem UTR Length per gene
#'
#' Calculate avergae length of tandem 3'UTRs 
#' @param object Seurat object
#' @param assay Name of polyAsite.Assay 
#' @param gene Name of gene with tandem UTR. 
#' @param pc Pseudocount to add for tandem UTR length calculation (default 0). 
#' @param N.cells Number of cells
#' @param Min.Reads.Per.Cell Use genes with minimum average read count (default 1).
#'
#' @return Returns a vector with average tandem UTR length
#' @export
#' @concept utilities

calcUTRlength <- function(gene,df,object,assay = assay, pc = 0, Min.Reads.Per.Cell = 1 ) {
  dsub <- subset(df, genes==gene)
  if ( sum(dsub$total)/ncol(object) > Min.Reads.Per.Cell){
    x=object[[assay]]@counts[rownames(dsub),]
    x1 = apply(x,2,FUN = function(z,PC=pc){(z+PC)/sum(z+PC)})
    if(all(rownames(x1) == rownames(dsub))){
      x2 = ((dsub$tandem_utr_reference_length * dsub$TandemUTRfraction) * (x1*100) ) 
      x3 = colSums(x2)/100
    }
    return(x3)
  }
}