
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
