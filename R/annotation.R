

#' Overlap PAS with polyA db 
#' @param object Seurat object containing a polyAsiteAssay
#' @param assay Name of polyAsiteAssay
#' @param polyAdb.file Location to find polyAdbv3 file 
#' @param max.dist Keep sites within this distance to a polyAdbv3 site. Default is 50 nucleotides. 
#' 
#' @importFrom GenomicRanges makeGRangesFromDataFrame 
#' @importFrom IRanges distanceToNearest trim
#' @importFrom plyranges anchor_3p mutate
#' @importFrom GenomeInfoDb seqlevelsStyle
#' 
#' @return Seurat object with annotated meta.features for polyAsiteAssay. Keeps all features 
#' in original matrix, but annotations columns will be populated with NAs if a feature does not 
#' overlap with a polyAdbv3 polyA site within the distance specified my max.dist.
#' @concept annotation
#' @export
#' 
GetPolyADbAnnotation <- function(
  object, 
  assay = "polyA",
  polyAdb.file = NULL, #
  max.dist = 50) 
{
  #readin polyAdbv3 and make GRanges file 
  if (!( file.exists(polyAdb.file))) {
    stop("Please check that you have specified the location of polyAdb.file correctly")
  }
  
  anno <- read.table(file = polyAdb.file, header = TRUE)
  GR.polyA.db = makeGRangesFromDataFrame( anno, 
                                          keep.extra.columns = TRUE, 
                                          seqnames.field = "hg38_Chromosome_format",
                                          start.field = "hg38_Position", 
                                          end.field = "hg38_Position",
                                          strand.field = "Strand")
  
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
  ranges = slot(object = object[[assay]], name = "ranges")
  if ( ! "strand" %in% colnames(object[[assay]]@meta.features)){
    object[[assay]] <- AddMetaData(object = object[[assay]] , metadata = as.character(strand(ranges)) , col.name = "strand")
  }
  #why include this? 
  #BiocGenerics::strand(ranges) <- S4Vectors::Rle(object[[assay]][["strand"]][,1])
  #mcols(ranges)$rn <- rownames(object[[assay]])
  if (  "*" %in% unique(strand(ranges)) ){
    if( table(strand(ranges))["*"] != 0 ){
      stop("Exiting. Cannot annotate unstranded PAS. Please remove unstranded PAS.")
    }
  } 
  if( !all(seqlevelsStyle(ranges) == seqlevelsStyle(GR.polyA.db)) ) {
    if( seqlevelsStyle(ranges)[1] == "UCSC"){
      warning("\n Annotation does not match between ranges and polyAdb.\n Annotation set to USCS.\n")
      seqlevelsStyle(GR.polyA.db) <- "UCSC"
    }else{
      warning("\n Annotation does not match between ranges and polyAdb.\n Annotation set to Ensembl.\n")
      seqlevelsStyle(GR.polyA.db) <- "Ensembl"
    }
  }
  GR = anchor_3p(ranges)
  GR_cleavage = mutate(GR, width = 1)
  #GR = trim(GR)
  #gr = ranges %>% anchor_3p() %>% mutate(width=1) %>% trim()
  OL = suppressWarnings(distanceToNearest(x = GR_cleavage, subject = GR.polyA.db,
                                          ignore.strand=FALSE))
  keep <- dplyr::filter(as.data.frame(OL), distance<= max.dist)
  
  #write this more to handle if htere are 2 sites equidistant
  if (sum(duplicated(queryHits(OL))) >0 ) {
    dups <- OL[duplicated(queryHits(OL))]
    
  }

  #re-write this 
  peak.df <- data.frame(GR)
  peak.df <- peak.df[,c("seqnames", "start", "end", "width", "strand")]
  peak.df$peak <- rownames(object[[assay]])

  #add information about what feature in object mathces polyAdbv3 peaks
  tmp <- peak.df[keep$queryHits,]
  tmp2 <- data.frame(GR.polyA.db)[keep$subjectHits,]
  tmp2 <- tmp2[,!(names(tmp2) %in%  c("strand", "seqnames", "start", "end", "width"))]
  tmp3 <- cbind(tmp, tmp2)
  
  meta.new <- left_join(peak.df, tmp3, by = c("seqnames", "start", "end", "strand", "width", "peak"))
  rownames(meta.new) <- meta.new$peak
  for( i in 1:ncol(meta.new)){
    object[[assay]] <- AddMetaData(object = object[[assay]] , metadata = meta.new[,colnames(meta.new)[i]] , col.name = colnames(meta.new)[i] )
  }
  return(object)
}
