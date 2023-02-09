
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Find polyA sites that are differentially used across cells, utilizing 
#' polyA residuals.
#' @param object An object
#' @param assay name of polyAsite assay to test. Default is polyA. 
#' @param ident.1 Identity class to find polyA sites for
#' @param ident.2 A second identity class for comparison. 
#' @param features Features to test. Default is all features.
#' @param covariates Vector of covariates to include in linear model.
#' @param gene.names Column name providing gene annotation of each polyA site.
#' Default is "Gene_Symbol"
#'
#' @importFrom stats lm relevel
#' 
#' @rdname FindDifferentialPolyA
#' @concept differential_polyA
#' @export
#'
FindDifferentialPolyA <- function(
    object,
    assay = "polyA",
    ident.1 = NULL,
    ident.2 = NULL,
    features = NULL,
    covariates = NULL,
    gene.names = "Gene_Symbol") {
  
  if (class(object[[assay]]) != "polyAsiteAssay") {
    stop("assay is not a polyAsite assay")
  }

  if (dim(GetAssayData(object, slot="scale.data", assay = assay))[1] == 0)  {
    stop ("No features found in scale.data slot for specified assay. Run CalcPolyAResiduals prior to FindPolyASites")
  }

  if (!(gene.names %in% colnames(object[[assay]][[]]))) {
    stop("Gene.names column not found in meta.features, please make sure
         you are specific gene.names correctly")
  }
  
  if (!(ident.1 %in% unique(Idents(object)))) {
    stop("ident.1 not found in object")
  }
  
  for (i in 1:length(covariates)) {
    if (is.na(match(covariates[[i]], colnames(object[[]])))) {
      stop("Covariates are not found in meta data, please make sure you are specifying correctly.")
    }
  }
  features <- features %||% rownames(x = object[[assay]]@scale.data)

  #check percentage difference
  percent.diff <- percentage.difference(object, cells.1 = WhichCells(object, idents = ident.1 ),
                                        cells.2 = WhichCells(object, idents = ident.2),
                                        features = features,
                                        assay=assay)

  df <- data.frame(idents = Idents(object))
  for (i in 1:length(covariates)) {
    df[,i+1] <- object[[]][,match(covariates[[i]], colnames(object[[]]))]
  }
  colnames(df) <- c("ident", covariates)
  
  r.matrix <- object[[assay]]@scale.data
  df$ident <- relevel(df$ident, ref = ident.2)     

  sub <- subset(df, ident %in% c(ident.1, ident.2))
  r.matrix.sub <- r.matrix[,rownames(sub)]

  all.models <- lapply(
    X = 1:nrow(x = r.matrix.sub),
    FUN = function(i) {
      sub$residuals <- as.numeric(r.matrix.sub[i, ])
      model <- lm(residuals ~  .,  data=sub)
      
      to_return <- data.frame(summary(model)$coefficients)
      to_return$coefficients <- rownames(to_return)
      colnames(to_return) <- c("Estimate", "std_error", "t", "p.value", "coefficient")
      to_return$peak <- rownames(r.matrix.sub)[i]
      return(to_return)
    }
  )
  results <- do.call(rbind, all.models)
  main.effects <- results[grep("ident", results$coefficient),]
  
  gene.idx = match(gene.names, colnames(object[[assay]][[]]))
  main.effects$symbol <- object[[assay]][[]][main.effects$peak,gene.idx]
  main.effects$percent.1 <- percentage.usage(object, 
                                             assay = assay, 
                                             cells = WhichCells(object, idents = ident.1), 
                                             features = main.effects$peak,
                                             gene.names = gene.names)
  main.effects$percent.2 <- percentage.usage(object, 
                                             assay = assay, 
                                             cells = WhichCells(object, idents = ident.2), 
                                             features = main.effects$peak,
                                             gene.names = gene.names)
  main.effects$p_val_adj <- main.effects$p.value * nrow(object[[assay]]@scale.data) 
  main.effects$p_val_adj[main.effects$p_val_adj  > 1] <- 1
  
  rownames(main.effects) <- main.effects$peak
  main.effects.return <- main.effects[,c("Estimate", "p.value", "p_val_adj", "percent.1", "percent.2", "symbol")]
  
  #order by p-value
  main.effects.return <- main.effects.return[ with(main.effects.return, order(p_val_adj, -Estimate)),]
  
  return(main.effects.return)
}


percentage.usage <- function( object,
                              assay = "polyA",
                              cells = NULL,
                              features = NULL,
                              gene.names = "Gene_Symbol") {
  
  df <- data.frame(peak=features)
  meta <- object[[assay]]@meta.features
  df$symbol <- meta[features,gene.names]
  
  genes.idx = match(gene.names, colnames(meta))
  genes.tmp <- meta[,genes.idx]
  peaks.tmp <- rownames(meta[genes.tmp %in% df$symbol,])

  df2 <- data.frame(peak=peaks.tmp)
  df2$counts1 <- rowSums(object[[assay]]@counts[peaks.tmp, cells])
  df2$gene <-meta[df2$peak, genes.idx]
  sum1 <- aggregate(df2$counts1, by=list(gene=df2$gene), FUN=sum)
  colnames(sum1) <- c("gene", "sum")
  df2 <- merge(df2, sum1, by="gene")
  df2$frac <- df2$counts1/df2$sum
  df2 <- df2[df2$peak %in% features,]
  return(as.vector(df2$frac))
}



percentage.difference <- function(object, cells.1, cells.2, features, assay) {
  df <- data.frame(peak=features)
  meta <- object[[assay]]@meta.features
  df$symbol <- meta[features,"Gene_Symbol"]
  
  peaks.tmp <- rownames(dplyr::filter(meta, Gene_Symbol %in% df$symbol))
  
  df2 <- data.frame(peak=peaks.tmp)
  df2$counts1 <- rowSums(object[[assay]]@counts[peaks.tmp, cells.1])
  df2$counts2 <- rowSums(object[[assay]]@counts[peaks.tmp, cells.2])
  df2$gene <-meta[df2$peak, "Gene_Symbol"]
  
  sum1 <- aggregate(df2$counts1, by=list(gene=df2$gene), FUN=sum)
  sum2 <- aggregate(df2$counts2, by=list(gene=df2$gene), FUN=sum)
  colnames(sum1)[2] <- "sum1"
  colnames(sum2)[2] <- "sum2"
  
  df2 <- merge(df2, sum1, by="gene")
  df2 <- merge(df2, sum2, by="gene")
  
  df2$frac1 <- df2$counts1/df2$sum1
  df2$frac2 <- df2$counts2/df2$sum2
  
  df2$percentage_diff <- df2$frac1 - df2$frac2
  rownames(df2) <- df2$peak
  df <- df2[features,"percentage_diff"]
  names(df) <- features
  return(as.vector(df))
}


