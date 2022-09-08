
#' Annotate PAS from GTF file 
#' @param object 
#' 
#' @importFrom plyr mapvalues
#' @importFrom GenomicRanges strand start end
#' @importFrom GenomeInfoDb renameSeqlevels genome seqlevelsStyle seqlevels
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom rtracklayer import
#' @export
#' 
#' @concept annotation

AnnotatePASfromGTF <- function (
    object, 
    assay,
    gtf.file,
    genome = genome, 
    hard_extension = 10, 
    extension = 2000,
    invert_strand = FALSE, 
    annotationType = "any", 
    transcriptDetails = TRUE,
    SequenceAnalysis = TRUE,
    PAS_Types = TRUE,
    pA_motif_max_position = 60, 
    AAA_motif_min_position = 10, 
    polystretch_length = 13, 
    max_mismatch = 1) 
{
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
  strand(ranges) <- Rle(object[[assay]][["strand"]][,1])
  mcols(ranges)$rn <- rownames(object[[assay]])
  if ("M" %in% seqlevels(ranges)) {
    ranges <- renameSeqlevels(x = ranges , value = plyr::mapvalues(seqlevels(ranges), from = "M", to = "MT"))
  }
  if ( table(strand(ranges))["*"] != 0){
    stop("Exiting. Cannot annotate unstranded PAS. Please remove unstranded PAS.")
  }
  gtf_gr <- import(gtf.file)
  gtf_TxDb <- makeTxDbFromGFF(gtf.file, format = "gtf")
  if( length(unique(genome(object[[assay]]))) == 0 ){
    stop("Exiting. No genome information provided in polyAsiteAssay.")
  }
  if( unique(genome(object[[assay]])) != unique(genome(genome))){
    stop("Exiting. Genome information provided in polyAsietAssay does not match BSgenome")
  }
  if( !all(seqlevelsStyle(gtf_TxDb) == seqlevelsStyle(genome)) | 
      !all(seqlevelsStyle(ranges) == seqlevelsStyle(genome)) |
      !all(seqlevelsStyle(ranges) == seqlevelsStyle(gtf_gr))) {
    
    if( seqlevelsStyle(ranges) == "UCSC"){
      warning("\nGenome annotation does not match between ranges, GTF and BSgenome.\n Annotation set to USCS.\n")
      seqlevelsStyle(gtf_TxDb) <- "UCSC"
      seqlevelsStyle(genome) <- "UCSC"
      seqlevelsStyle(gtf_gr) <- "UCSC"
    }else{
      warning("\nGenome annotation does not match between ranges, GTF and BSgenome.\n Annotation set to Ensembl.\n")
      seqlevelsStyle(ranges) <- "Ensembl"
      seqlevelsStyle(gtf_TxDb) <- "Ensembl"
      seqlevelsStyle(genome) <- "Ensembl"
      seqlevelsStyle(gtf_gr) <- "Ensembl"
    }
  }
  
  # ranges <- ranges[which(as.character(seqnames(ranges)) == "X")[1:50]]
  
  # Modified from Sierra package
  annot.df <- annotate_gr_from_gtf(gr = ranges, gtf_gr = gtf_gr, 
                                   gtf_TxDb = gtf_TxDb, genome = genome, hard_extension = hard_extension, extension = extension,  
                                   invert_strand = invert_strand, annotationType = annotationType, transcriptDetails = transcriptDetails, 
                                   SequenceAnalysis = SequenceAnalysis, PAS_Types = PAS_Types,
                                   pA_motif_max_position = pA_motif_max_position, AAA_motif_min_position = AAA_motif_min_position,
                                   polystretch_length = polystretch_length, max_mismatch = max_mismatch)
  # annot.df$rn = mcols(ranges)$rn
  
  expected_cols = c("symbol", "gene_id" , "gene_biotype", "UTR3", "UTR3_extension", 
                    "UTR5", "intron", "exon", "CDS","pA_motifs","pA_motif_pos",
                    "pA_stretch_len","pA_stretch_pos","misprime",
                    "cleavage_coord" ,"cds_end" ,"tandem_utr_start" ,
                    "tandem_utr_end" ,"tandem_utr_reference_length",
                    "distance2UTRstart","TandemUTRfraction","is.tandem" )
  idx = which(colnames(annot.df) %in% expected_cols)
  if( length(idx) > 0){
    
    # matching is switched off
    # due to the possibility of peaks that are different on strand but share the same coordinates, this currently breaks
    # Since all entries remain in the same order, I for now cbind columns without checking matching  names
    #m = match( rownames(object[[assay]]@meta.features), annot.df$rn)
    #if( all( rownames(object[[assay]]@meta.features) == annot.df$rn[m]) ){
    #  for( i in idx){
    #    object[[assay]] <- AddMetaData(object = object[[assay]] , metadata = annot.df[,colnames(annot.df)[i]][m] , col.name = colnames(annot.df)[i] )
    #  }
    #}
    #else{
    #  stop("Exiting. Feature names do not match.")
    #}
    
    for( i in idx){
      object[[assay]] <- AddMetaData(object = object[[assay]] , metadata = annot.df[,colnames(annot.df)[i]] , col.name = colnames(annot.df)[i] )
    }
    
    
  }
  else{
    stop("\nExiting. No annotation returned.\n")
  }
  
  return(object)
  
}




#' Annotate Granges from GTF
#' 
#' importFrom
# This function has been modified from Sierra package
# The modification allows for single nucleotide mapping
# of PAS assuming that the PAS 3'end indicates the true
# transcript end. Moreover, this function now allows
# for mapping to extended 3'UTRs if no other gene is
# annotated. And, we now automatically classify PAS
# to reside as tandem PAS or alternative last exons.

# TO DO: 1) shift gtf_txdb by including hard and soft 
# extensions to allow for mapping of PAS downstream 
# of annotated UTRs. 2) Annotate tandem PAS and ALEs.
# 3) Distance to stop codon. 4) APA motif detail.


annotate_gr_from_gtf = function(gr, gtf_gr = NULL, gtf_TxDb, genome = NULL, 
                                hard_extension = 10, extension = 2000, invert_strand = FALSE,
                                annotationType = "any", transcriptDetails = TRUE,
                                SequenceAnalysis = TRUE, PAS_Types = FALSE,
                                pA_motif_max_position = 60, AAA_motif_min_position = 10, 
                                polystretch_length = 10, max_mismatch = 1) 
{
  if (is.null(gtf_gr)) {
    warning("No gtf file provided.\n")
    return(NULL)
  }
  if (invert_strand) {
    gr <- invertStrand(gr)
  }
  GenomicRanges::mcols(gtf_gr) <- GenomicRanges::mcols(gtf_gr)[c("type", "gene_id", "gene_name","gene_biotype")]
  if (length(intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(gtf_gr))) == 0) {
    warning("Exiting. Genome anotation information does not match.\nseqlevelsStyle set to \"Ensembl\".")
    seqlevelsStyle(gr) <- "Ensembl"
    seqlevelsStyle(gtf_gr) <- "Ensembl"
  }
  if (length(grep(pattern = "type", x = colnames(as.data.frame(gtf_gr)))) == 0) {
    message("Annotation reference does not have a type field. This field is used\nto extract gene, transcript, exon, UTR information. Cannot continue...\n")
    return(NULL)
  }
  listed_annotations <- names(table(gtf_gr$type))
  if (length(grep(pattern = "gene", x = listed_annotations)) > 0) {
    genes_gr <- gtf_gr[gtf_gr$type == "gene"]
  }
  else {
    message("Reference has no recognised labels to annotate. Cannot continue... \n")
    return(NULL)
  }
  # map extensions to genes
  seq_ranges <- seqinfo(genome) %>% as("GRanges")
  seq_stranded_ranges <- c(
    plyranges::mutate(seq_ranges, strand = "+"),
    plyranges::mutate(seq_ranges, strand = "-"))
  # get gaps not covered by gene features
  gaps_gr <- suppressWarnings(plyranges::setdiff_ranges_directed(seq_stranded_ranges, genes_gr))
  
  # get transcript extensions
  extensions_gr <- gaps_gr %>%
    plyranges::join_overlap_inner_directed(flank_downstream(genes_gr,1)) %>%
    plyranges::anchor_5p() %>%
    plyranges::mutate(width=pmin(width, extension))
  
  backup.gr = gr
  # mutate ranges to single nt anchored on 3'end to improve annotation of transcript 3'ends
  # for hard extension, we assume that the PAS is shifted slightly downstream of its true location
  # instead of extending all reference ranges by the value of hard_extension, 
  # we shift the end coordinate of the query gr by hard_extension towards its 5'end
  if(hard_extension != 0){
    gr = gr %>% plyranges::anchor_3p() %>% plyranges::mutate(width=1) %>% shift_upstream(shift = hard_extension) %>% IRanges::trim()
  } else{
    gr = gr %>% plyranges::anchor_3p() %>% plyranges::mutate(width=1) %>% IRanges::trim()
  }
  
  annotate_info <- gene_Labels(gr = gr, reference_gr = genes_gr, annotationType = annotationType)
  extension_info <- gene_Labels(gr = gr, reference_gr = extensions_gr, annotationType = annotationType)
  if (is.null(annotate_info)) {
    warning("tbd")
    return(NULL)
  }
  
  df <- as.data.frame(backup.gr)
  df$symbol <- ""
  df$symbol[annotate_info$idx_to_annotate] <- annotate_info$identified_gene_symbols
  if (extension > 0){
    df$symbol[extension_info$idx_to_annotate] <- extension_info$identified_gene_symbols
  }
  
  # add ensemble id and biotype
  m = match(df$symbol , mcols(genes_gr)$gene_name)
  df$gene_id <- mcols(genes_gr)$gene_id[m]
  df$gene_biotype <- mcols(genes_gr)$gene_biotype[m]
  
  
  if (transcriptDetails) {
    
    genes_gr <- genes_gr %>% plyranges::mutate(gene_strand = strand)
    
    db_trans <- transcriptsBy(gtf_TxDb) %>%
      delist("gene_id") %>%
      dplyr::select(gene_id, tx_id, tx_name) %>%
      left_join_mcols(genes_gr, "gene_id") %>%
      dplyr::filter(strand == gene_strand) %>% #Forbid antisense isoforms
      dplyr::select(tx_id, tx_name, gene_id, gene_name, gene_biotype)
    
    
    cat("\nAnnotating 3'UTRs")
    df$UTR3 <- ""
    UTR_3_GR <- GenomicFeatures::threeUTRsByTranscript(gtf_TxDb, use.names = TRUE) %>%
      delist("tx_name") %>%
      plyranges::select(tx_name) %>%
      left_join_mcols(db_trans, "tx_name")
    all_UTR_3_hits <- GenomicAlignments::findOverlaps(gr,UTR_3_GR, type = annotationType)
    idx_to_annotate_3UTR <- S4Vectors::queryHits(all_UTR_3_hits)
    df$UTR3[idx_to_annotate_3UTR] <- "YES"
    
    if (extension > 0){
      cat("\nAnnotating 3'UTR extensions")
      extensions_gr <- gaps_gr %>%
        plyranges::join_overlap_inner_directed(flank_downstream(UTR_3_GR,1)) %>%
        plyranges::anchor_5p() %>%
        plyranges::mutate(width=pmin(width, extension))
      all_UTR3_extension_hits <- GenomicAlignments::findOverlaps(gr,extensions_gr, type = annotationType)
      idx_to_annotate_3UTRextension <- S4Vectors::queryHits(all_UTR3_extension_hits)
      df$UTR3_extension <- ""
      df$UTR3_extension[idx_to_annotate_3UTRextension] <- "YES"
    }
    
    
    cat("\nAnnotating 5'UTRs")
    df$UTR5 <- ""
    UTR_5_GR <- GenomicFeatures::fiveUTRsByTranscript(gtf_TxDb)
    all_UTR_5_hits <- GenomicAlignments::findOverlaps(gr, UTR_5_GR, type = annotationType)
    identified_5UTRs <- UTR_5_GR[S4Vectors::subjectHits(all_UTR_5_hits)]$exon_id
    idx_to_annotate_5UTR <- S4Vectors::queryHits(all_UTR_5_hits)
    
    df$UTR5[idx_to_annotate_5UTR] <- "YES"
    cat("\nAnnotating introns")
    df$intron <- ""
    introns_GR <- GenomicFeatures::intronsByTranscript(gtf_TxDb)
    all_intron_hits <- GenomicAlignments::findOverlaps(gr, introns_GR, type = annotationType)
    identified_introns <- introns_GR[S4Vectors::subjectHits(all_intron_hits)]
    idx_to_annotate_introns <- S4Vectors::queryHits(all_intron_hits)
    
    df$intron[idx_to_annotate_introns] <- "YES"
    cat("\nAnnotating exons")
    my_exons <- gtf_gr[gtf_gr$type == "exon"]
    all_exon_hits <- GenomicAlignments::findOverlaps(gr, my_exons, type = annotationType)
    idx_to_annotate_exons <- S4Vectors::queryHits(all_exon_hits)
    df$exon <- ""
    
    df$exon[idx_to_annotate_exons] <- "YES"
    cat("\nAnnotating CDS")
    my_CDS <- gtf_gr[gtf_gr$type == "CDS"]
    all_CDS_hits <- GenomicAlignments::findOverlaps(gr,  my_CDS, type = annotationType)
    idx_to_annotate_CDS <- S4Vectors::queryHits(all_CDS_hits)
    df$CDS <- ""
    df$CDS[idx_to_annotate_CDS] <- "YES"
  }
  
  
  if(SequenceAnalysis){
    
    if (!is.null(genome)) {
      cat("\nAnalysing sequence motifs around cleavage sites\n")
      if (isS4(genome)) {
        
        # batch mode, faster then iterating one by one
        motif_details <- BaseCompositionFast(genome, coord = backup.gr, mismatch = max_mismatch, 
                                             PAS_uptream_length = pA_motif_max_position, 
                                             offset = 0, 
                                             PAS_downstream_length = 25,
                                             AT_max = 20, AT_min = 5)
        
        df$pA_motifs <- motif_details$pA_motif
        df$pA_motif_pos <- motif_details$pA_motif_pos
        df$pA_stretch_len <- motif_details$pA_stretch_len
        df$pA_stretch_pos <- motif_details$pA_stretch_pos
        df$misprime <- ifelse( motif_details$pA_stretch_pos <= AAA_motif_min_position & motif_details$pA_stretch_len >= polystretch_length , TRUE, FALSE )
        df$misprime[is.na(motif_details$pA_stretch_pos)==TRUE] <- FALSE
        
        # iterating one by one, slower, depricated
        #pb <- txtProgressBar(min = 0, max = length(backup.gr), style = 3)
        #motif_details <- lapply(X = 1:length(backup.gr), FUN = function(i,backup.gr, pb) {
        #  setTxtProgressBar(pb, i)
        #  nextgr <- backup.gr[i]
        #  BaseComposition(genome, coord = nextgr, mismatch = max_mismatch, 
        #                  PAS_uptream_length = pA_motif_max_position)
        #}, backup.gr, pb)
        # cat("\n")
        # df$pA_motifs <- sapply( motif_details, FUN = function(x) { 
        #   x$pA_motif
        # })
        # df$pA_motif_pos <- sapply( motif_details, FUN = function(x) { 
        #   x$pA_motif_pos
        # })
        # df$pA_stretch_len <- sapply( motif_details, FUN = function(x) { 
        #   x$pA_stretch_len
        # })
        # df$pA_stretch_pos <- sapply( motif_details, FUN = function(x) { 
        #   x$pA_stretch_pos
        # })
        # df$misprime <- sapply( motif_details, FUN = function(x) { 
        #   if( is.na( x$pA_stretch_pos ) ){
        #     FALSE
        #   }else{
        #     if ( x$pA_stretch_pos <= AAA_motif_min_position & x$pA_stretch_len >= polystretch_length ) {
        #       TRUE
        #     }else{
        #       FALSE
        #     }
        #   }
        # })
      }  
      else {
        warning("Genome object is not a BSgenome S4 object.\n
              Cannot annotate for sequence motifs")
      }
    }
  }
  
  if(PAS_Types){
    
    cat("\nCreate tandem 3'UTR reference from GTF (~10min)")
    tandem.gr = MakeTandemUTRReferenceFromGTF(txdb = gtf_TxDb, genes_gr = genes_gr, gaps_gr = gaps_gr , extension = extension )
    
    # annotate PAS by overlapping
    cat("\nAnnotating tandem 3'UTRs\n")
    tandem_hits <- GenomicAlignments::findOverlaps(gr,tandem.gr)
    q.idx = S4Vectors::queryHits(tandem_hits)
    s.idx = S4Vectors::subjectHits(tandem_hits)
    
    gr$cds_end = NA
    gr$tandem_utr_start = NA
    gr$tandem_utr_end = NA
    gr$tandem_utr_reference_length = NA
    gr$distance2UTRstart = NA
    gr$TandemUTRfraction = NA
    gr$cleavage_coord = NA
    
    gr$cds_end[q.idx] <- tandem.gr$cds_end[s.idx]
    gr$tandem_utr_start[q.idx] <- anchor_5p(tandem.gr[s.idx]) %>% plyranges::mutate(width=1) %>% start()
    gr$tandem_utr_end[q.idx] <- anchor_3p(tandem.gr[s.idx]) %>% plyranges::mutate(width=1) %>% end()
    gr$tandem_utr_reference_length <- abs(gr$tandem_utr_end - gr$tandem_utr_start)
    gr$distance2UTRstart <- abs(start(gr) - gr$tandem_utr_start)
    gr$TandemUTRfraction <- gr$distance2UTRstart / gr$tandem_utr_reference_length
    gr$cleavage_coord <- end(anchor_3p(gr))
    # find tandem UTRs with > 1 cleavage site overlap
    d = s.idx[which( duplicated(s.idx) == T )]
    tandem.q.idx = q.idx[which(s.idx %in% d)]
    gr$is.tandem = FALSE
    gr$is.tandem[tandem.q.idx] = TRUE
    
    
    #match back to output df
    m = match(df$rownames , gr$rownames)
    df$cleavage_coord <- gr$cleavage_coord[m]
    df$cds_end <- gr$cds_end[m]
    df$tandem_utr_start <- gr$tandem_utr_start[m]
    df$tandem_utr_end <- gr$tandem_utr_end[m]
    df$tandem_utr_reference_length <- gr$tandem_utr_reference_length[m]
    df$distance2UTRstart <- gr$distance2UTRstart[m]
    df$TandemUTRfraction <- gr$TandemUTRfraction[m]
    df$is.tandem <- gr$is.tandem[m]
    
  }
  
  return(df)
}


## helper functions
hard_extend <- function(GR)
  GR %>% anchor_5p() %>% plyranges::mutate(width=width+hard_extension) %>% trim()

just_mcols <- function(GR)
  GR %>% mcols() %>% as.data.frame()

left_join_mcols <- function(left, right, by) {
  df <- just_mcols(left) %>%
    plyranges::mutate(.row. = row_number()) %>%
    dplyr::left_join(just_mcols(right), by=by)
  
  result <- left[df$.row.,]
  mcols(result) <- dplyr::select(df, -.row.)
  result
}

delist <- function(grl, id_col) {
  GR <- unlist(grl)
  mcols(GR)[[id_col]] <- names(GR)
  names(GR) <- NULL
  GR
}


MakeTandemUTRReferenceFromGTF = function(txdb = gtf_TxDb, genes_gr = NULL, gaps_gr = NULL , extension = 0 ){
  
  genes_gr <- genes_gr %>% plyranges::mutate(gene_strand = strand)
  
  db_trans <- transcriptsBy(txdb) %>%
    delist("gene_id") %>%
    dplyr::select(gene_id, tx_id, tx_name) %>%
    left_join_mcols(genes_gr, "gene_id") %>%
    dplyr::filter(strand == gene_strand) %>% #Forbid antisense isoforms
    dplyr::select(tx_id, tx_name, gene_id, gene_name, gene_biotype)
  
  # stop codon (stop codon is the last 3 nt in cds)
  cds_by_tx <- cdsBy(txdb, by="tx", use.names=TRUE) %>% delist("tx_name")
  plus = split(cds_by_tx, strand(cds_by_tx))[["+"]]
  minus = split(cds_by_tx, strand(cds_by_tx))[["-"]]
  cds.plus = split(plus, plus$tx_name) 
  cds.minus = split(minus, minus$tx_name) 
  cds.ends.plus = sapply( end(cds.plus) , max )
  cds.ends.minus = sapply( start(cds.minus) , min )
  cds.ends = c(cds.ends.plus,cds.ends.minus)
  
  #Make 3'UTR db
  Tx = GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE)
  UTR_3_GR <- Tx %>%
    delist("tx_name") %>%
    plyranges::select(tx_name) %>%
    left_join_mcols(db_trans, "tx_name")
  UTR_3_GR$UID <- c(1:length(UTR_3_GR))
  
  if (extension > 0){
    extensions_gr <- gaps_gr %>%
      plyranges::join_overlap_inner_directed(flank_downstream(UTR_3_GR,1)) %>%
      plyranges::anchor_5p() %>%
      plyranges::mutate(width=pmin(width, extension))
    idx = which( UTR_3_GR$UID %in% extensions_gr$UID)
    ext = UTR_3_GR[idx]
    if (!all(ext$UID == extensions_gr$UID)){
      stop("exiting. 3'UTR UIDs don't match.")
    }
    ext = ext %>% anchor_5p() %>% plyranges::mutate(width=width+width(extensions_gr)) %>% trim()
    UTR_3_GR = UTR_3_GR[-idx]
    UTR_3_GR = c(UTR_3_GR,ext)
    UTR_3_GR = UTR_3_GR[  order(UTR_3_GR$UID, decreasing = F) ]
  }
  
  # remove genes without unique strand
  plus = UTR_3_GR[which(strand(UTR_3_GR) == "+")] 
  minus = UTR_3_GR[which(strand(UTR_3_GR) == "-")] 
  if ( length( unique(plus[which(plus$gene_name %in% minus$gene_name)]$gene_name)) > 0){
    g = unique(plus[which(plus$gene_name %in% minus$gene_name)]$gene_name)
    UTR_3_GR = UTR_3_GR[-which(UTR_3_GR$gene_name %in% g)]
  }
  
  #split by strand
  UTR_3_GR.plus = split(UTR_3_GR, strand(UTR_3_GR))[["+"]]
  UTR_3_GR.minus = split(UTR_3_GR, strand(UTR_3_GR))[["-"]]
  
  
  #1
  # a) find most proximal 3'UTR positions for each transcript, 
  # b) most distal 3' UTR position for each gene,
  # c) and lengths of all annotated 3'UTRs
  # this includes hard extension if applicable
  starts.plus = UTR_3_GR.plus %>% group_by( tx_name ) %>% summarise(utr_start = min(start))
  ends.plus = UTR_3_GR.plus %>% group_by( tx_name ) %>% summarise(utr_end = max(end))
  starts.minus = UTR_3_GR.minus %>% group_by( tx_name ) %>% summarise(utr_start = max(end))
  ends.minus = UTR_3_GR.minus %>% group_by( tx_name ) %>% summarise(utr_end = min(start))
  starts = rbind(starts.plus,starts.minus)
  ends = rbind(ends.plus,ends.minus)
  # width
  Tx = split(UTR_3_GR, UTR_3_GR$tx_name) #makes sure it includes any 3'UTR extensions of terminal exons
  W = sum(width(Tx)) 
  # exons
  ExonLengths = sapply( width(Tx) , FUN = function(z){paste0( z, collapse = ';' )})
  # Get 3'UTR introns
  Tx.plus = Tx[UTR_3_GR.plus$tx_name]
  Tx.minus = Tx[UTR_3_GR.minus$tx_name]
  Tx.plus.mult = Tx.plus[which( elementNROWS(Tx.plus) > 1)]
  Tx.minus.mult = Tx.minus[which( elementNROWS(Tx.minus) > 1)]
  
  Int.plus = lapply( Tx.plus.mult , FUN = function(x){ gaps( x , start = min(start(x))) })
  Int.plus = as(Int.plus,"CompressedGRangesList" )
  Int.width.plus = sapply( width(Int.plus) , FUN = function(z){paste0( z, collapse = ';' )})
  Ex.starts.plus = sapply( start(Tx.plus) , FUN = function(z){paste0( z, collapse = ';' )})
  Ex.ends.plus = sapply( end(Tx.plus) , FUN = function(z){paste0( z, collapse = ';' )})
  
  Int.minus = lapply( Tx.minus.mult , FUN = function(x){ gaps( x , start = min(start(x))) })
  Int.minus = as(Int.minus,"CompressedGRangesList" )
  Int.width.minus = sapply( width(Int.minus) , FUN = function(z){paste0( rev(z), collapse = ';' )})
  Ex.starts.minus = sapply( end(Tx.minus) , FUN = function(z){paste0( z, collapse = ';' )})
  Ex.ends.minus = sapply( start(Tx.minus) , FUN = function(z){paste0( z, collapse = ';' )})
  
  Int.width = c(Int.width.plus,Int.width.minus)
  Ex.starts = c(Ex.starts.plus,Ex.starts.minus)
  Ex.ends = c(Ex.ends.plus,Ex.ends.minus)
  
  # complete UTR reference
  DF = unique(as.data.frame(UTR_3_GR)[c("gene_name","gene_id","tx_name","seqnames","strand")])
  m = match(DF$tx_name, names(cds.ends))
  DF$cds_end = cds.ends[m]
  m = match(DF$tx_name,starts$tx_name)
  DF$utr_start = starts$utr_start[m]
  m = match(DF$tx_name,ends$tx_name)
  DF$utr_end = ends$utr_end[m]
  m = match(DF$tx_name,names(W))
  DF$utr_length = W[m]
  m = match(DF$tx_name,names(ExonLengths))
  DF$utr_exon_width = ExonLengths[m]
  m = match(DF$tx_name,names(Int.width))
  DF$utr_intron_width = Int.width[m]
  m = match(DF$tx_name,names(Ex.starts))
  DF$utr_exon_starts = Ex.starts[m]
  m = match(DF$tx_name,names(Ex.ends))
  DF$utr_exon_end = Ex.ends[m]
  
  # a) find most proximal 3'UTR positions for each transcript
  # done (=utr_start)
  # b) most distal 3' UTR position for each gene
  G = split(DF, f = DF$gene_name)
  gene_utr_ends = sapply(G , FUN = function(z){
    if(nrow(z) == 1){
      z$utr_end
    }else{
      if(unique(z$strand) == "+"){
        max(z$utr_end)
      }else{
        min(z$utr_end)
      }
    }
  })
  m = match(DF$gene_name, names(gene_utr_ends))
  DF$gene_utr_end = gene_utr_ends[m]
  # c) and lengths of all annotated 3'UTRs
  # done (= utr_length)
  
  #2
  # choose representative (R) transcript with longest UTR for each unique 3' utr start position (S)
  # split by 3'utr start
  # select most distal 3' end
  # all TX with same coord will be stored in tx_name
  S = split(DF, f = paste0(DF$gene_name,"_",DF$utr_start))
  R = lapply(S , FUN = function(z){
    if(nrow(z) == 1){
      z
    }else{
      if(unique(as.character(z$strand)) == "+"){
        idx = which(z$utr_end == max(z$utr_end))
      }else{
        idx = which(z$utr_end == min(z$utr_end))
      }
      if(length(idx) > 1){
        z$tx_name[idx] <- paste0( z$tx_name[idx], collapse = ";")
      }
      z[idx[1],]
    }
  })
  DF.representative = do.call(rbind.data.frame, R)
  rownames(DF.representative) <- NULL
  
  #3
  #find position of most distal stop codon ( = 3'utr start minus 1) per gene, and corresponding ID for representative transcript harboring it
  #compile a list of representative transcript IDs
  # split by gene
  G = split(DF.representative, f = DF.representative$gene_name)
  R = lapply(G , FUN = function(z){
    if(nrow(z) == 1){
      z$representative_utr_end <- z$utr_end
      z$representative_utr_start <- z$utr_start
      z
    }else{
      
      if(unique(as.character(z$strand)) == "+"){
        idx.u = which(z$utr_end == max(z$utr_end))[1]
        z$representative_utr_end <- z$utr_end[idx.u]
        idx = which(z$utr_start == max(z$utr_start))
        z$representative_utr_start <- z$utr_start[idx[1]]
      }else{
        idx.u = which(z$utr_end == min(z$utr_end))[1]
        z$representative_utr_end <- z$utr_end[idx.u]
        idx = which(z$utr_start == min(z$utr_start))
        z$representative_utr_start <- z$utr_start[idx[1]]
      }
      z
    }
  })
  DF.representative = do.call(rbind.data.frame, R)
  rownames(DF.representative) <- NULL
  
  #4 select representative UTR
  # select maximal UTR (greatest dist between start and end)
  G = split(DF.representative , f = DF.representative$gene_name)
  R = lapply(G , FUN = function(z){
    if(nrow(z) == 1){
      z
    }else{
      if(unique(as.character(z$strand)) == "+"){
        
        l = z$representative_utr_end - z$representative_utr_start
        idx = which(l == max(l))
        if(length(idx)>1){
          #use most distal CDS end
          z = z[idx,]
          idx = which(z$cds_end == max(z$cds_end))[1]
          z[idx,]
        }else{
          z[idx,]
        }
        
      }else{
        l = z$representative_utr_start - z$representative_utr_end
        idx = which(l == max(l))
        if(length(idx)>1){
          #use most distal CDS end
          z = z[idx,]
          idx = which(z$cds_end == max(z$cds_end))[1]
          z[idx,]
        }else{
          z[idx,]
        }
      }
    }
  })
  
  DF.representative = do.call(rbind.data.frame, R)
  rownames(DF.representative) <- NULL
  
  
  plus = makeGRangesFromDataFrame(subset(DF.representative, strand == "+"),
                                  keep.extra.columns=TRUE,
                                  ignore.strand=FALSE,
                                  seqnames.field=c("seqnames", "seqname",
                                                   "chromosome", "chrom",
                                                   "chr", "chromosome_name",
                                                   "seqid"),
                                  start.field="representative_utr_start",
                                  end.field=c("representative_utr_end"),
                                  strand.field="strand",
                                  starts.in.df.are.0based=FALSE)
  minus = makeGRangesFromDataFrame(subset(DF.representative, strand == "-"),
                                   keep.extra.columns=TRUE,
                                   ignore.strand=FALSE,
                                   seqnames.field=c("seqnames", "seqname",
                                                    "chromosome", "chrom",
                                                    "chr", "chromosome_name",
                                                    "seqid"),
                                   start.field="representative_utr_end",
                                   end.field=c("representative_utr_start"),
                                   strand.field="strand",
                                   starts.in.df.are.0based=FALSE)
  tandem.gr = c(plus,minus)
  
  return(tandem.gr)
  
}


################################################################################################
## Function taken from Sierra
## This function has been designed to be called from annotate_gr_from_gtf
gene_Labels <- function(gr, reference_gr, annotationType)
{
  if ( "M" %in% seqlevels(gr) & "MT" %in% seqlevels(reference_gr) ) {
    gr <- renameSeqlevels(x = gr , value = plyr::mapvalues(seqlevels(gr), from = "M", to = "MT"))
  }
  if ( "MT" %in% seqlevels(gr) & "M" %in% seqlevels(reference_gr) ) {
    gr <- renameSeqlevels(x = gr , value = plyr::mapvalues(seqlevels(gr), from = "MT", to = "M"))
  }
  
  all_hits <- GenomicAlignments::findOverlaps(gr , reference_gr, type = annotationType)
  if (length(reference_gr) == 0)
  { warning("No entries in reference to annotate. Cannot continue")
    return(NULL)
  }
  if (length(all_hits) == 0)
  { sanityCheck <- length(intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(reference_gr)))
  msg <- paste0("No peaks aligned to any entry within gtf reference.",
                "\n Sanity check: ", sanityCheck,
                " seqnames (i.e. chromosomes) match between peak and reference file")
  warning(msg)
  return(NULL)
  }
  
  identified_gene_symbols <- reference_gr[S4Vectors::subjectHits(all_hits)]$gene_name
  idx_to_annotate <- S4Vectors::queryHits(all_hits)
  
  multi_annotations <- which(table(idx_to_annotate) > 1)
  unique_annotations <- unique(idx_to_annotate)
  
  multi_gene_IDs <- unlist(sapply(names(multi_annotations),FUN=function(x) {
    newID <- paste(unique(identified_gene_symbols[which(x== idx_to_annotate)]),collapse=",")
    # rep(newID, length(which(x== idx_to_annotate)))
  }))
  multi_idx <- as.numeric(names(multi_gene_IDs))  # These are the indexes to annotate
  
  to_convert <- lapply(multi_idx,FUN = function(x) {which(idx_to_annotate == x)})
  
  if (length(to_convert) > 0)
  {   for(i in 1:length(to_convert))
  {
    identified_gene_symbols[to_convert[[i]]] <- multi_gene_IDs[[i]]
  }
  }
  
  return(list (identified_gene_symbols=identified_gene_symbols, idx_to_annotate=idx_to_annotate))
}



#############################################################################
# Function inspired by Sierra package, but simplified
# We assume the 3' coordinate marks the 3'end of the transcript

BaseComposition <- function (genome = NULL, coord = NULL, offset = 0, 
                             PAS_downstream_length = 25,
                             PAS_uptream_length = 60, 
                             mismatch = 1, AT_max = 20, AT_min = 5) 
{
  if (!isS4(genome)) {
    warning("genome is not a BSgenome S4 object. Cannot continue.")
    return(NULL)
  }
  if( !class(coord) == "GRanges" ){
    stop("Exiting. Expecting object of class GRanges for \"coord\".")
  }
  
  gr.PAS = coord %>% plyranges::anchor_3p() %>% plyranges::mutate(width=1) %>% IRanges::trim()
  gr.upstream = gr.PAS %>% plyranges::anchor_3p() %>% plyranges::mutate(width=PAS_uptream_length)
  gr.downstream = gr.PAS %>% plyranges::shift_downstream(shift = 1) %>% plyranges::anchor_5p() %>% plyranges::mutate(width=PAS_downstream_length)
  
  seq.upstream <- BSgenome::getSeq(genome, gr.upstream )
  seq.downstream <- BSgenome::getSeq(genome, gr.downstream )
  
  # check for cleavage motifs
  PAS_signals = c("AATAAA", "ATTAAA", "TATAAA", "AGTAAA", "AATACA", 
                  "CATAAA", "AATATA", "GATAAA", "AATGAA", "AATAAT", "AAGAAA", 
                  "ACTAAA", "AATAGA", "ATTACA", "AACAAA", "ATTATA", "AACAAG", 
                  "AATAAG")
  signal_matches = lapply( PAS_signals , FUN = function( pattern , 
                                                         seq = seq.upstream){
    Biostrings::vmatchPattern(pattern = pattern, 
                              subject = seq)[[1]]
  })
  names(signal_matches) = PAS_signals
  idx = which( sapply(signal_matches , length ) > 0)
  if(length(idx)>0){
    # motifs (s) and match position (p) are in matching order 
    # separated by semicolon. 
    # Match position is relative to 3'end (=PAS).
    signal_matches = signal_matches[idx]
    s = paste0( rep( names(signal_matches), times = sapply(signal_matches,length)), collapse = ";" )
    p = sapply(signal_matches , FUN = function(j, len = PAS_uptream_length){ 
      len - start(j) + 1 }) %>% unlist() %>% paste0( collapse = ";")
  }else{
    s = NA
    p = NA 
  }
  
  # check for mis priming (= polyA-stretch downstream of PAS)
  pA = Biostrings::DNAString(paste0(rep("A",AT_max),collapse = ""))
  PA_signals = as.character( Views(subject = pA , start = 1, end = seq(AT_max,AT_min,-1)))
  
  signal_matches = list()
  for (i in seq_along(PA_signals)){
    signal_matches[[i]] <- Biostrings::vmatchPattern(pattern = PA_signals[i], subject = seq.downstream, 
                                                     max.mismatch = mismatch)[[1]]
    names(signal_matches)[i] <- PA_signals[i]
    if( length(signal_matches[[i]]) > 0){
      break
    }
  }
  
  idx = which( sapply(signal_matches , length ) > 0)
  if(length(idx)>0){
    # selecting the longest polyA stretch
    # returning its length in seq.downstream
    signal_matches = signal_matches[idx[1]]
    l = width(signal_matches[[1]])[1]
    m = start(signal_matches[[1]])[1]
  }else{
    m=NA
    l=0
  }
  
  return(list(pA_motif = s,
              pA_motif_pos = p,
              pA_stretch_len = l, 
              pA_stretch_pos = m,
              sequence = list(upstream = seq.upstream, downstream = seq.downstream)))
  
}


BaseCompositionFast <- function (genome = NULL, coord = NULL, offset = 0, 
                                 PAS_downstream_length = 25,
                                 PAS_uptream_length = 60, 
                                 mismatch = 1, AT_max = 20, AT_min = 5) 
{
  if (!isS4(genome)) {
    warning("genome is not a BSgenome S4 object. Cannot continue.")
    return(NULL)
  }
  if( !class(coord) == "GRanges" ){
    stop("Exiting. Expecting object of class GRanges for \"coord\".")
  }
  
  gr.PAS = coord %>% plyranges::anchor_3p() %>% plyranges::mutate(width=1) %>% IRanges::trim()
  gr.upstream = gr.PAS %>% plyranges::anchor_3p() %>% plyranges::mutate(width=PAS_uptream_length)
  gr.downstream = gr.PAS %>% plyranges::shift_downstream(shift = 1) %>% plyranges::anchor_5p() %>% plyranges::mutate(width=PAS_downstream_length)
  
  seq.upstream <- BSgenome::getSeq(genome, gr.upstream )
  seq.downstream <- BSgenome::getSeq(genome, gr.downstream )
  
  # check for cleavage motifs
  PAS_signals = c("AATAAA", "ATTAAA", "TATAAA", "AGTAAA", "AATACA", 
                  "CATAAA", "AATATA", "GATAAA", "AATGAA", "AATAAT", "AAGAAA", 
                  "ACTAAA", "AATAGA", "ATTACA", "AACAAA", "ATTATA", "AACAAG", 
                  "AATAAG")
  
  
  signal_matches = lapply( PAS_signals , FUN = function( pattern, seq = seq.upstream){
    Biostrings::vmatchPattern(pattern = pattern, subject = seq)
  })
  names(signal_matches) = PAS_signals
  batchmatches = list()
  for(i in 1:length(signal_matches)){
    batchmatches[[i]] = GetBatchMatches(SM = signal_matches[[i]], cse = names(signal_matches)[i] , len = PAS_uptream_length )
  }
  names(batchmatches) = names(signal_matches)
  
  M = do.call(cbind.data.frame,lapply(batchmatches,FUN = function(j){j[[1]]}))
  s = apply(M,1,FUN = function(j){gsub(":+",";",gsub("^:+|:+$","",paste0(j,collapse = ":")))})
  M = do.call(cbind.data.frame,lapply(batchmatches,FUN = function(j){j[[2]]}))
  p = apply(M,1,FUN = function(j){gsub(":+",";",gsub("^:+|:+$","",paste0(j,collapse = ":")))})
  
  
  
  # check for mis priming (= polyA-stretch downstream of PAS)
  pA = Biostrings::DNAString(paste0(rep("A",AT_max),collapse = ""))
  PA_signals = as.character( Views(subject = pA , start = 1, end = seq(AT_max,AT_min,-1)))
  
  signal_matches = lapply( PA_signals , FUN = function( pattern, seq = seq.downstream){
    Biostrings::vmatchPattern(pattern = pattern, subject = seq,  max.mismatch = mismatch)
  })
  names(signal_matches) = PA_signals
  batchmatches = list()
  for(i in 1:length(signal_matches)){
    batchmatches[[i]] = GetBatchMatchesPA(SM = signal_matches[[i]], pa = names(signal_matches)[i])
  }
  names(batchmatches) = names(signal_matches)
  
  
  M = do.call(cbind.data.frame,lapply(batchmatches,FUN = function(j){j[[1]]}))
  l = apply(M,1,FUN = function(j){ nchar( strsplit( gsub("^;+","",gsub("NA","", gsub(":+",";",gsub("^:+|:+$","",paste0(j,collapse = ":"))))) , ";")[[1]][1]) })
  M = do.call(cbind.data.frame,lapply(batchmatches,FUN = function(j){j[[2]]}))
  m = apply(M,1,FUN = function(j){  strsplit( gsub("^;+","",gsub("NA","", gsub(":+",";",gsub("^:+|:+$","",paste0(j,collapse = ":"))))) , ";")[[1]][1] })
  
  
  
  
  
  return(list(pA_motif = s,
              pA_motif_pos = p,
              pA_stretch_len = l, 
              pA_stretch_pos = m,
              sequence = list(upstream = seq.upstream, downstream = seq.downstream)))
  
}

GetBatchMatches = function(SM = signal_matches[[1]], cse = names(signal_matches)[1] , len = PAS_uptream_length ){
  idx = which( sapply(SM , length ) > 0)
  l = sapply(SM , length )
  
  if(length(idx)>0){
    # motifs (s) and match position (p) are in matching order 
    # separated by semicolon. 
    # Match position is relative to 3'end (=PAS).
    s = sapply( SM, FUN = function(sm,CSE = cse){paste0(rep(CSE,times = length(sm)) , collapse = ";")})
    p = sapply( SM, FUN = function(sm,LEN = len){ (len - start(sm) + 1)  %>% paste0( collapse = ";") } )
    
  }else{
    s = NA
    p = NA 
  }
  return(list(s,p))
}

GetBatchMatchesPA = function(SM = signal_matches[[1]], pa = names(signal_matches)[1] ){
  idx = which( sapply(SM , length ) > 0)
  l = sapply(SM , length )
  
  if(length(idx)>0){
    # motifs (s) and match position (p) are in matching order 
    # separated by semicolon. 
    # Match position is relative to 3'end (=PAS).
    s = sapply( SM, FUN = function(sm,PA = pa){paste0(rep(PA,times = length(sm)) , collapse = ";")})
    p = sapply( SM, FUN = function(sm){ start(sm)  %>% paste0( collapse = ";") } )
    
  }else{
    s = NA
    p = NA 
  }
  return(list(s,p))
}

##############################

# This function classifies PAS in UTRs as single PAS, tandem PAS or 
# alternative last exon (ALE) PAS.
# This is a separate function, as we assume users may filter their polyAsiteAssay
# based on basic annotation. Filtering will have an impact on PAS classification 
# in the types described above. 

# For this analysis, we group 3'UTRs by start coordinate and reduce the ranges
# into the longest continuous 3'UTR (plus extension). ALEs are defined by 
# alternative 3'UTR five-prime end. if a PAS overlaps multiple ALEs, we assign 
# the most three-prime. This allows us that each PAS is assigned to only 
# one transcript/3'UTR isoform. We remove all transcripts that are not annotated 
# as protein coding for all protein coding genes (e.g. non-sense mediated decay).

# In addition, we return the distance of each PAS/3'UTR end to 3'UTR start 
# (=3'UTR length). 

# How to deal with spliced 3'UTRs

GetUTRtypes <- function(
    object, 
    assay,
    gtf.file,
    genome = NULL,
    chrLen_file = NULL,
    hard_extension = 10, 
    extension = 2000
){
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
  strand(ranges) <- Rle(object[[assay]][["strand"]][,1])
  mcols(ranges)$rn <- rownames(object[[assay]])
  if ( table(strand(ranges))["*"] != 0){
    stop("Exiting. Cannot annotate unstranded PAS. Please remove unstranded PAS.")
  }
  gtf_gr <- rtracklayer::import(gtf.file)
  gtf_TxDb <- GenomicFeatures::makeTxDbFromGFF(gtf.file, format = "gtf")
  if( length(unique(genome(object[[assay]]))) == 0 ){
    stop("Exiting. No genome information provided in polyAsiteAssay.")
  }
  if ("M" %in% seqlevels(ranges)) {
    ranges <- renameSeqlevels(x = ranges , value = plyr::mapvalues(seqlevels(ranges), from = "M", to = "MT"))
  }
  if (is.null(genome) & is.null(chrLen_file)){
    stop("Exiting. Please provide a BSgenome or chromosome length file")
  }
  if(!is.null(genome)){
    if( unique(genome(object[[assay]])) != unique(GenomeInfoDb::genome(genome))){
      stop("Exiting. Genome information provided in polyAsietAssay does not match BSgenome")
    }
    seq_ranges <- seqinfo(genome) %>% as("GRanges")
  }else{
    warning("ChrLen_file does not provide genome information. Proceeding without genome verification.")
    ChrLen = read.delim(chrLen_file , header = FALSE, sep = "\t", stringsAsFactors = F )
    seq_ranges <-  GRanges(seqnames = as.character(ChrLen$V1),
                           ranges = IRanges(rep(1,nrow(ChrLen)), end = ChrLen$V2, names = as.character(ChrLen$V1)),
                           strand = Rle(strand("*"), nrow(ChrLen) ))
  }
  if( !all(seqlevelsStyle(gtf_TxDb) == seqlevelsStyle(seq_ranges)) &
      !all(seqlevelsStyle(ranges) == seqlevelsStyle(seq_ranges))) {
    warning("Genome annotation does not match between ranges, GTF and BSgenome.\nGenome annotation set to Ensembl.")
    seqlevelsStyle(ranges) <- "Ensembl"
    seqlevelsStyle(gtf_TxDb) <- "Ensembl"
    seqlevelsStyle(seq_ranges) <- "Ensembl"
  }
  
  
  ### map extensions to genes
  gr = ranges 
  # get annotated gene ranges
  GenomicRanges::mcols(gtf_gr) <- GenomicRanges::mcols(gtf_gr)[c("type", "gene_id", "gene_name","gene_biotype","transcript_id","transcript_biotype")]
  if (length(grep(pattern = "type", x = colnames(as.data.frame(gtf_gr)))) == 0) {
    message("Annotation reference does not have a type field.\nThis field is used\nto extract gene, transcript, exon, UTR information. Cannot continue")
    return(NULL)
  }
  
  backup.gr = gr # ranges backup before transforming to single nt resolution
  # mutate ranges to single nt anchored on 3'end to improve annotation of transcript 3'ends
  # for hard extension, we assume that the PAS is shifted slightly downstream of its true location
  # instead of extending all reference ranges by the value of hard_extension, 
  # we shift the end coordinate of the query gr by hard_extension towards its five-prime end
  if(hard_extension != 0){
    gr = gr %>% plyranges::anchor_3p() %>% plyranges::mutate(width=1) %>% plyranges::shift_upstream(shift = hard_extension) %>% IRanges::trim()
  } else{
    gr = gr %>% plyranges::anchor_3p() %>% plyranges::mutate(width=1) %>% IRanges::trim()
  }
  
  
  # step 1
  # build a transcript database for efficient filtering and tandem-PA identification
  # keep all protein-coding transcripts
  # add info for: 
  # stop codon coord, 3'UTR start, 3'UTR end (incl. pot. extension), 
  # 3'UTR length, semicolon-separated exon starts, semicolon-separated exon ends
  
  # Build 3'UTR database
  # functions
  hard_extend <- function(GR)
    GR %>% anchor_5p() %>% plyranges::mutate(width=width+hard_extension) %>% trim()
  
  just_mcols <- function(GR)
    GR %>% mcols() %>% as.data.frame()
  
  left_join_mcols <- function(left, right, by) {
    df <- just_mcols(left) %>%
      plyranges::mutate(.row. = row_number()) %>%
      dplyr::left_join(just_mcols(right), by=by)
    
    result <- left[df$.row.,]
    mcols(result) <- dplyr::select(df, -.row.)
    result
  }
  
  delist <- function(grl, id_col) {
    GR <- unlist(grl)
    mcols(GR)[[id_col]] <- names(GR)
    names(GR) <- NULL
    GR
  }
  
  genes_gr <- gtf_gr[gtf_gr$type == "gene"] %>%
    plyranges::mutate(gene_strand = strand)
  
  tx_gr <- gtf_gr[gtf_gr$type == "transcript"] %>%
    plyranges::mutate(gene_strand = strand)
  # For protein coding genes, remove transcripts that are not protein coding
  tx_gr <- tx_gr[which(tx_gr$gene_biotype == "protein_coding" & tx_gr$transcript_biotype == "protein_coding")]
  
  # get chromosome ranges
  seq_stranded_ranges <- c(
    plyranges::mutate(seq_ranges, strand = "+"),
    plyranges::mutate(seq_ranges, strand = "-"))
  
  # get transcripts
  db_trans <- GenomicFeatures::transcriptsBy(gtf_TxDb) %>%
    delist("gene_id") %>%
    dplyr::select(gene_id, tx_id, tx_name) %>%
    left_join_mcols(genes_gr, "gene_id") %>%
    dplyr::filter(strand == gene_strand) %>% #Forbid antisense isoforms
    dplyr::filter( tx_name %in% tx_gr$transcript_id) %>% #Allow only protein coding tx for protein coding genes
    dplyr::select(tx_id, tx_name, gene_id, gene_name, gene_biotype)
  
  
  # get 3'UTRs
  utr3s <- GenomicFeatures::threeUTRsByTranscript(gtf_TxDb, use.names = TRUE)
  db_utr3s <- utr3s %>%
    delist("tx_name") %>%
    plyranges::select(tx_name) %>%
    dplyr::filter( tx_name %in% tx_gr$transcript_id) %>% #Allow only protein coding tx for protein coding genes
    left_join_mcols(db_trans, "tx_name")
  
  # gaps == what is not a transcript
  db_gaps <- plyranges::setdiff_ranges_directed(seq_stranded_ranges, db_trans)
  
  # make transcript/3'utr extension
  db_ext <- db_gaps %>%
    plyranges::join_overlap_inner_directed(plyranges::flank_downstream(db_utr3s,1)) %>%
    plyranges::anchor_5p() %>%
    plyranges::mutate(width=pmin(width, extension)) 
  db_ext_uniq = unique(db_ext)
  
  # extent 3'UTRs by extension
  OL = findOverlaps( query = db_utr3s , subject = flank_upstream(db_ext_uniq,1))
  db_utr3s_ext = db_utr3s
  db_utr3s_ext[queryHits(OL)] = db_utr3s[queryHits(OL)] %>% anchor_5p() %>% plyranges::mutate(width= width + width(db_ext_uniq[subjectHits(OL)]) ) 
  
  # add if exon is adjacent to stop codon
  # get stop codons
  stop_gr <- gtf_gr[gtf_gr$type == "stop_codon"] %>%
    plyranges::mutate(gene_strand = strand) %>%
    plyranges::filter(transcript_id %in% db_utr3s_ext$tx_name) %>%
    plyranges::filter( width > 2)
  
  OL = findOverlaps( query = flank_upstream(db_utr3s_ext,1) , subject = stop_gr)
  db_utr3s_ext$stop.adjacent <- FALSE
  db_utr3s_ext$stop.adjacent[unique(queryHits(OL))] <- TRUE
  
  stop_gr = stop_gr %>% plyranges::anchor_5p() %>% plyranges::mutate(width=1) %>% IRanges::trim()
  
  # get 3'UTR start, first filter for 5'exon by start coord
  db_utr3s_ext_5primeEnds = db_utr3s_ext %>% plyranges::anchor_5p() %>% plyranges::mutate(width=1) %>% plyranges::mutate(s=start)
  db_utr3s_ext_filt = db_utr3s_ext_5primeEnds %>% 
    group_by(tx_name) %>% 
    plyranges::filter( s == min(s)) %>% #this doesn't work for both strands
    ungroup() 
  
  # get 3'UTR ends, first filter for e'exon by end coord
  db_utr3s_ext_3primeEnds = db_utr3s_ext %>% plyranges::anchor_3p() %>% plyranges::mutate(width=1) %>% plyranges::mutate(s=start)
  db_utr3s_ext_filt = db_utr3s_ext_3primeEnds %>% 
    group_by(tx_name) %>% 
    plyranges::filter( s == max(s)) %>%  #this doesn't work for both strands
    ungroup() 
  
  
  #1
  #find most a) proximal 3'UTR positions for each transcript, 
  # b) most distal 3' UTR position for each gene, 
  # and c) lengths of all annotated 3'UTRs
  
  # add to db_trans: stop coding coord, 3'UTR start, 3'UTR end
  db_trans$stop_codon = NA
  db_trans$utr3_start = NA
  db_trans$utr3_end = NA
  m = match(db_trans$tx_name , stop_gr$transcript_id)
  db_trans$stop_codon <- start(stop_gr)[m]
  m = match(db_trans$tx_name , db_utr3s_ext_filt$tx_name)
  db_trans$utr3_start <- start(db_utr3s_ext_filt)[m]
  
  
  
  
  
  # add additional info for filtering
  db_utr3s_ext$utr.start.coord  <- start(db_utr3s_ext %>% plyranges::anchor_start() %>% plyranges::mutate(width=1))
  db_utr3s_ext = db_utr3s_ext %>% plyranges::mutate(utr.length=width)
  # choose representative transcript with longest UTR for each unique 3'utr start position
  db_utr3s_ext_filt = db_utr3s_ext %>% 
    group_by(utr.start.coord) %>%
    plyranges::filter( utr.length == max(utr.length)) %>%
    ungroup() %>%
    unique()
  
  
  G = split(db_utr3s_ext_filt, db_utr3s_ext_filt$tx_id)
  
  # find overlaps, given the PAS
  OL = findOverlaps( query = gr , subject = db_utr3s_ext_filt)
  
  # find duplicated subjectHits, indicative for >1 pA-site overlapping the same 3'UTR fragment ( = tandem UTR)
  dups = subjectHits(OL)[ which(duplicated(subjectHits(OL))) ]
  dups.subjectHits.idx = which( subjectHits(OL) %in% dups)
  gr$is.tandem = FALSE
  gr$is.tandem[queryHits(OL)[dups.subjectHits.idx]] = TRUE
  
}
