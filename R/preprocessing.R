
#' The function reads in output from polyApipe and returns a count matrix
#' @param counts.file file containing the counts from polyApipe in .tsv format
#' @param peaks.file gff file that requires the seqid, start, end, and strand of each peak
#' @param sep Separators
#' @param filter.chromosomes If TRUE (default), only include main chromosomes. 
#' @param min.features Cell must have at least this many feature to be included.
#' @param min.cells A feature must be present in at least this manhy cells.
#' @param verbose Print output.
#'
#' @return Returns count matrix with peaks in the form "chromosome:start,end:strand"
#'  such as "10:100000560,100000859,100000560,100000859:+"
#'
#'
#' @importFrom data.table fread
#' @importFrom Matrix sparseMatrix colSums rowSums
#' @importFrom rtracklayer readGFF
#'
#' @export
#' @concept preprocessing
ReadPolyApipe <- function(counts.file, peaks.file = NULL, sep = c(":",",",":"),
                          filter.chromosomes = TRUE, min.features = NULL, min.cells = NULL , 
                          verbose=TRUE)  {

  if (verbose) {
    message("Reading count file ...")
  }
  the_counts <- fread( file = counts.file, sep = "\t", header=TRUE,
                                   showProgress = verbose,
                                   colClasses = c("factor","factor", "numeric"),
                                   nrows=-1)
  colnames(the_counts) <- c("peak","cell","count")

  # First construct the sparse matrix of all the counts
  if (verbose) {
    message("Creating sparse matrix ...")
  }
  counts_matrix <- sparseMatrix(
    i = as.integer(the_counts$peak),
    j = as.integer(the_counts$cell),
    x = as.integer(the_counts$count))
  colnames(counts_matrix) <- levels(the_counts$cell)
  rownames(counts_matrix) <- levels(the_counts$peak)

  if(! is.null( peaks.file )){
    if (verbose) {
      message("Reading peak file ...")
    }
    the_peaks = readGFF( filepath = peaks.file)
    the_peaks$peak <- gsub( '"', '', the_peaks$peak)
    m = match( rownames(counts_matrix)  , the_peaks$peak )
    if ( length (which( is.na(m) == TRUE )) > 0 ){
      if ( length(which( is.na(m) == TRUE )) == nrow(counts_matrix) ){
        stop("Exiting. peaks.file does not match count_matrix rownames.")
      }else{
        n = nrow(counts_matrix)
        del = which(! rownames(counts_matrix) %in% the_peaks$peak )
        idx = which( rownames(counts_matrix) %in% the_peaks$peak )
        counts_matrix = counts_matrix[idx,]
        warning(paste0( "Removed ",length(del)," of ", n ," peaks not present in peaks.file"))
        m = match( rownames(counts_matrix)  , the_peaks$peak ) #update

        the_peaks$newRN = apply( the_peaks[,c("seqid","start", "end", "strand")], 1,FUN =
                                   function(x){ paste0(x[1], sep[1],x[2], sep[2], x[3], sep[3], x[4]   )} )


      }
    }else{
      the_peaks$newRN = apply( the_peaks[,c("seqid","start", "end", "strand")], 1,FUN =
                                 function(x){ paste0(x[1], sep[1],x[2], sep[2], x[3], sep[3], x[4]    )} )
      rownames(counts_matrix) <- the_peaks$newRN[m]
    }
  }
  
  if(!is.null(min.features)){
    rs = rowSums( counts_matrix > 0 )
    del = which(rs < min.features)
    keep = which(rs >= min.features)
    counts_matrix = counts_matrix[keep,]
    if (verbose) {
      message(paste0("Removed ",length(del)," feature(s) covered in less than ",min.features," cells ..."))
    }
  }
  
  if(!is.null(min.cells)){
    cs = colSums( counts_matrix > 0 )
    del = which(cs < min.cells)
    keep = which(cs >= min.cells)
    counts_matrix = counts_matrix[,keep]
    if (verbose) {
      message(paste0("Removed ",length(del)," cell(s) with less than ",min.cells," features ..."))
    }
  }
  
  if (filter.chromosomes) {
    chr = sapply(rownames(counts_matrix) , FUN = function(j) { strsplit(j,":")[[1]][1]})
    accept = c(paste0("^",c(c(1:23),"X","Y","MT","M")),
               paste0("chr",c(c(1:23),"X","Y","MT","M")))
    del = grep(paste0(accept,collapse = "|"),chr, invert = T)
    keep = grep(paste0(accept,collapse = "|"),chr, invert = F)
    counts_matrix = counts_matrix[keep,]
    if (verbose) {
      message("Removed ", length(del), " entries falling on scaffolds")
    }
  }
  
  if (verbose) {
  message("Loaded ", nrow(counts_matrix), " x ", ncol(counts_matrix), " feature by cell matrix")
  }
  
  return(counts_matrix)
}



#' Finds variable features based on polyA residuals
#' @param object Object containing a polyAsite assay
#' @param nfeatures Number of features to select as top variable features.
#' @param gene.names Column in meta features that contains gene names. At most one feature per gene will be selected.
#' @param selection.method How to calculate polyA residuals. If "residuals" (default), 
#' will rank all polyA sites by their variance and pick at most 1 polyA site per gene. 
#' Otherwise, will use Seurat FindVariableFeatures functions.
#' 
#' @rdname FindVariableFeatures
#' @importFrom utils head tail
#' @concept preprocessing
#' @export
#' @method FindVariableFeatures polyAsiteAssay
#' @seealso \code{\link[Seurat]{FindVariableFeatures}} for additional parameters.

FindVariableFeatures.polyAsiteAssay <- function(
    object,
    nfeatures = 2000,
    selection.method = "residuals",
    gene.names = "symbol",
    ...)
  {
  if (selection.method == "residuals") {
    
    if (dim(GetAssayData(object, slot="scale.data"))[1] == 0)  {
      stop ("No features found in scale.data slot. Run CalcPolyAResiduals prior to FindVariableFeatures.")
    }
    
    if (!(gene.names %in% colnames(object[[]]))) {
      stop("Gene.names column not found in meta.features, please make sure 
         you are specific gene.names correctly")
    }
    
    var <- data.frame(var = apply(object@scale.data, 1, function(x) var(x[x != 0])))
    var$symbol <- object[[]][rownames(var), gene.names]
    var <- var[order(var$var, decreasing = TRUE), ] #sort by maximum variance 
    var.unique <- var[!duplicated(var$symbol),]
    
    nfeatures <- min(nfeatures, nrow(x = nrow(var.unique)))
  
    VariableFeatures(object) <- head(rownames(var.unique), n=nfeatures)
  } else {
    assay <- as(object = object[[DefaultAssay(object)]], Class = "ChromatinAssay")
    tmp <- Seurat::FindVariableFeatures(assay, nfeatures = nfeatures, selection.method = selection.method, ...)
    VariableFeatures(object) <- VariableFeatures(tmp)
  }
  return(object)
}



