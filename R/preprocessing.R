

#peaks_file = "~/nygc_cluster_hwessels/Science2020_Arunalacham/PAS_polyA_peaks.gff"

ReadPolyApipe <- function(counts_file, peaks_file = NULL, sep = c(":",",",":"), verbose=TRUE){
  #' The function reads in output from polyApipe and returns a count matrix
  #' @param counts_file file containing the counts from polyApipe in .tsv format
  #' @param peaks_file gff file that requires the seqid, start, end, and strand of each peak
  #' @param sep Separators
  #'
  #' @return Returns count matrix with peaks in the form "chromosome:start,end:strand"
  #'  such as "10:100000560,100000859,100000560,100000859:+"
  #'
  #'
  #' @importFrom data.table fread
  #' @importFrom Matrix sparseMatrix
  #' @importFrom rtracklayer readGFF
  #'
  #' @export
  #' @concept preprocessing

  # Not using read_tsv because of compression.(segfault.)
  # from polyApipe https://github.com/MonashBioinformaticsPlatform/polyApipe/blob/master/polyApiper/R/data_loading.R
  # NB: https://github.com/tidyverse/readr/issues/610
  if (verbose) {
    message("Reading count file ...")
  }
  the_counts <- fread( file = counts_file, sep = "\t", header=TRUE,
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

  if(! is.null( peaks_file )){
    if (verbose) {
      message("Reading peak file ...")
    }
    the_peaks = readGFF( filepath = peaks_file)
    the_peaks$peak <- gsub( '"', '', the_peaks$peak)
    m = match( rownames(counts_matrix)  , the_peaks$peak )
    if ( length (which( is.na(m) == TRUE )) > 0 ){
      if ( length(which( is.na(m) == TRUE )) == nrow(counts_matrix) ){
        stop("Exiting. peaks_file does not match count_matrix rownames.")
      }else{
        n = nrow(counts_matrix)
        del = which(! rownames(counts_matrix) %in% the_peaks$peak )
        idx = which( rownames(counts_matrix) %in% the_peaks$peak )
        counts_matrix = counts_matrix[idx,]
        warning(paste0( "Removed ",length(del)," of ", n ," peaks not present in peak_file"))
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

  message("Loaded ", nrow(counts_matrix), " x ", ncol(counts_matrix), " matrix of counts")
  return(counts_matrix)
}



#' Finds variable features based on polyA residuals
#' @param object Object containing a polyAsite assay
#' @param nfeatures Number of features to select as top variable features.

#' @rdname FindVariableFeatures
#'
#' @concept preprocessing
#' @export
#' @method FindVariableFeatures polyAsiteAssay

FindVariableFeatures.polyAsiteAssay <- function(
    object,
    nfeatures = 2000,
    ...
) {

  nfeatures <- min(nfeatures, nrow(x = nrow(object)))
  var <- apply(object@scale.data, 1, function(x) var(x[x != 0]))
  VariableFeatures(object) <- head(names(sort(var, decreasing = TRUE)), n=nfeatures)
  return(object)
}



