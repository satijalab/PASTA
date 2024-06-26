
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Calculate PolyA Residuals
#'
#'
#' @param object Seurat object containing a polyAsiteAssay
#' @param assay Name of polyAsiteAssay to be used in calculating polyAresiduals
#' @param features Features to include in calculation of polyA residuals.
#' Default is to use all features.
#' @param background Identity of cells to use as background.
#' Default is to use all cells as a background.
#' @param gene.names Column containing the gene where each polyA site is annotated.
#' Default is symbol.
#' @param min.counts.background Features with at least this many counts in the background cells are included in calculation
#' @param min.variance Sets minimum variance. Default is 0.1.
#' @param sample.n Max number of observations to sample in each bin when performing regularization. Default is 1000.
#' @param do.center Return the centered residuals. Default is TRUE.
#' @param do.scale Return the scaled residuals. Default is TRUE.
#' @param residuals.max Clip residuals above this value. Default is NULL (no clipping).
#' @param residuals.min Clip residuals below this value. Default is NULL (no clipping).
#' @param number.bins Number of bins used to perform regularization. Default is 30.
#' @param verbose Print messages.
#'
#'
#' @return Returns a Seurat object with polyAresiduals assay
#'
#' @export
#' @concept residuals
#'
CalcPolyAResiduals <- function(object,
                               assay="polyA",
                               features=NULL,
                               background = NULL,
                               gene.names = "Gene_Symbol",
                               min.counts.background = 5,
                               min.variance = 0.1,
                               sample.n = 1000,
                               do.scale = FALSE,
                               do.center = FALSE,
                               residuals.max = NULL,
                               residuals.min = NULL,
                               number.bins = 30,
                               verbose=TRUE)
 {
  if(verbose) {
    message("Calculating background distribution")
  }

  #if features in NULL, then specify all features in polyA assay
  if (is.null(features)) {
    features <- rownames(LayerData(object, assay=assay, layer="counts"))
  }

  #if background is NULl, then make a dummy variable
  if (is.null(background)) {
    message("Using all cells in order to estimate background distribution")
    object$dummy <- "all"
    Idents(object) <- object$dummy
    background.use = "all"
  }  else {
    background.use = background
    if (!(background %in% unique(Idents(object)))) {
      stop("background must be one of the Idents of seurat object")
    }
    message(paste0("Using ", background, " as background distribution"))
  }


  #check if symbols are contained in meta features
  if (!(gene.names %in% colnames(object[[assay]]@meta.features))) {
    stop("Gene.names column not found in meta.features, please make sure
         you are specific gene.names correctly")
  }

  if (sum(is.na(object[[assay]]@meta.features[features,gene.names])) > 0) {
    features.no.anno <- features[is.na(object[[assay]]@meta.features[features,gene.names])]
    message(paste0("Removing ", length(features.no.anno), " sites without a gene annotation"))
    features <- setdiff(features, features.no.anno)
  }

  ##############################################################################
  #get pseudobulked fraction of reads from background
  background.dist <- GetBackgroundDist(object = object, features = features,
                                       background = background.use,
                                       gene.names = gene.names,
                                       assay = assay,
                                       min.counts.background = min.counts.background)
  #remove features without gene annotation
  background.dist <- subset(background.dist, gene!="_-")
  background.dist <- subset(background.dist, gene!="_+")
  features.use <- background.dist$peak

  ##############################################################################
  #calculate sum of counts within each gene
  m <- LayerData(object = object, assay=assay, layer="counts")
  m <- m[background.dist$peak,]
  #m <- m[order(match(rownames(m), background.dist$peak)), ]
  gene.sum <- rowsum(m, group=background.dist$gene)
  genes <- rownames(gene.sum)

  ##############################################################################
  #fit dirichlet multinomial distribution
  if(verbose) {
    message("Running Dirichlet Multionmial Regression")
  }

  ncells = dim(object)[2]
  background.cells <- WhichCells(object, idents=background)
  m.background <- as.matrix(m[,background.cells], nrow = nrow(m))

  #fit dirichlet multionmial for each gene
  res <- lapply(genes, DirichletMultionmial, background.dist=background.dist,
                m.background = m.background, gene.sum=gene.sum,  ncells = ncells)

  ec <- data.frame(dummy = 1:length(colnames(object)))
  var <- data.frame(dummy = 1:length(colnames(object)))
  for(i in 1:length(res)) {
    if (!(is.null(res[[i]]))) {
      ec <- cbind(ec, res[[i]]$ec)
      var <- cbind(var, res[[i]]$var)
    }
  }

  ec <- t(ec[, -1])
  var <- t(var[, -1])
  colnames(ec) <- colnames(object)
  colnames(var) <- colnames(object)

  m <- m[rownames(ec),]

  ##############################################################################
  ### regularize dirichlet multinomial variance
  if(verbose) {
    message("Regularizing Dirichlet Multionmial Variance")
  }

  var.reg <- RegDMVar(ec = ec, var = var, m = m, m.background = m.background,
                      background.dist = background.dist,
                      gene.sum = gene.sum, background.cells = background.cells,
                      min.variance = min.variance,
                      number.bins = number.bins,
                      sample.n = sample.n)
  #calculate residual matrix
  residual.matrix <- (m-ec) / sqrt(var.reg)
  residual.matrix <- as.matrix(residual.matrix, nrow = nrow(residual.matrix))
  #M1 <- as(residual.matrix, "dgCMatrix")

  #change to same order as counts slot
  features.order <- rownames(residual.matrix)[order(
    match(
      rownames(residual.matrix),
      rownames(LayerData(object = object, assay=assay, layer="counts"))
    )
  )]
  residual.matrix <- residual.matrix[features.order,]
  residual.matrix <- scale(residual.matrix, center=do.center, scale= do.scale )
  if (!is.null(residuals.max)) {
    residual.matrix[residual.matrix > residuals.max] <- residuals.max
  }

  if (!is.null(residuals.min)) {
    residual.matrix[residual.matrix < residuals.min] <- residuals.min
  }
  #change default assay
  DefaultAssay(object = object) <- assay
  #need to met SetAssayData, GetAssayData for residuals
  LayerData(object, layer = "scale.data") <- residual.matrix
  object <- LogSeuratCommand(object = object)
}


#' Get Background Distribution
#'
#' Calculated Pseudobulk Ratios of Each Isoform within a gene for background distribution.
#' @param object Seurat object containing a polyAsiteAssay
#' @param assay Name of polyAsiteAssay to be used in calculating polyAresiduals
#' @param features Features to include in calculation of polyA residuals.
#' If NULL, use all features.
#' @param background Identity of cells to use as background
#' If NULL, uses all cells combined as a background
#' @param gene.names Name of column containing gene annotations
#' @param min.counts.background Features with at least this many counts in the background cells are included in calculation
#'
#' @return Returns a data frame containing all peaks within genes that have multiple polyA sites that meet min.counts.background criteria
#'
#' @importFrom stats aggregate
#' @concept residuals
#'
GetBackgroundDist <- function(object, features, background, gene.names, assay,  min.counts.background) {
  # returns the pseudobulked background distribution for peaks specified
  # must contain gene information in meta data

  suppressMessages(nt.pseudo <- AverageExpression(object, features = features, assays = assay, slot="counts"))
  nt.pseudo <- data.frame(background = nt.pseudo[[1]][,background]) #subset just the background
  nt.pseudo$background <- nt.pseudo$background * sum(Idents(object)==background)
  nt.pseudo$gene <- paste0(object[[assay]]@meta.features[features, gene.names], "_", object[[assay]]@meta.features[features, "strand"])
  nt.pseudo$peak <- rownames(nt.pseudo)

  nt.pseudo <- nt.pseudo[nt.pseudo$background>min.counts.background,] #subset to peaks with min number of counts
  genes.use <- nt.pseudo$gene[duplicated(nt.pseudo$gene)] #only use genes with at least 2 peaks per gene
  nt.pseudo <- nt.pseudo[nt.pseudo$gene %in% genes.use,]

  if ( length(genes.use)  ==  0) {
    stop("Found no genes with more than 2 features within a gene. Please make sure you are including
         all peaks within a gene you would like to include.")
  }

  tmp <- aggregate(nt.pseudo$background, list(nt.pseudo$gene), FUN=sum)
  colnames(tmp) <- c("gene", "sum")
  nt.pseudo <- merge(nt.pseudo, tmp, by="gene")
  nt.pseudo$frac <- nt.pseudo$background/ nt.pseudo$sum
  return(nt.pseudo)
}


#' Run Dirichlet Multinomial Distribution
#'
#' Fit dirichlet multinomial distribution on each peak within a gene using background cells.
#' Then calculate expected value and variance for each cell based on estimates from dirichlet multionimial regression.
#'
#' @param gene.test which gene to use
#' @param background.dist dataframe containing the isoform ratios for each
#' @param m.background matrix of background distribution
#' @param gene.sum sum of count within each gene for each cell
#' @param ncells number of cells
#'
#' @return Returns a list where first element is matrix of expected values for each peak within the genes,
#' second value is matrix of variance for each peak within the gene
#'
#' @importFrom MGLM MGLMfit
#' @concept residuals
#'
DirichletMultionmial <- function(
  gene.test,
  background.dist,
  m.background,
  gene.sum,
  ncells
) {
  peaks <- background.dist$peak[background.dist$gene == gene.test]
  t <- m.background[rownames(m.background) %in% peaks,]
  t <- t(t)
  fit <- try(compareFit <- suppressWarnings(MGLMfit(t, dist="DM")), silent=TRUE)

  if (!inherits(fit, "try-error")) {
    param <- compareFit@estimate
    sum.p <- sum(param)
    n <-  as.numeric(gene.sum[gene.test,])
    expect.tmp <- data.frame(matrix(nrow=ncells, ncol=length(peaks)))
    var.tmp <- data.frame(matrix(nrow=ncells, ncol=length(peaks)))
    #calculate expected and variance for each peak
    for (i in 1:length(peaks)) {
      expect.x <- n*param[i]/sum.p
      var.x <- expect.x*(1- param[i]/sum.p)*(n + sum.p)/(1+sum.p)
      expect.tmp[,i] <- expect.x
      var.tmp[,i] <- var.x
    }
    colnames(expect.tmp) <- peaks
    colnames(var.tmp) <- peaks
    return(list(ec = expect.tmp, var = var.tmp))
  }
}


#' Run Dirichlet Multionmial Distribution
#'
#' Calculated Pseudobulk Ratios of Each Isoform within a gene for background distribution.
#'
#'
#'
#' @return Returns a data frame containing all peaks within genes that have multiple polyA sites that meet min.counts.background criteria
#'
#' @importFrom dplyr left_join
#' @importFrom stats quantile
#' @importFrom gplm kreg
#' @concept residuals
#'
#'
RegDMVar <- function(ec,
                     var, m,
                     m.background,
                     background.dist,
                     gene.sum,
                     background.cells,
                     min.variance = min.variance,
                     number.bins = number.bins,
                     sample.n = 1000
                     ) {
  expected.counts.df <- data.frame(expected.counts = matrix(ec, ncol=1))
  expected.counts.df$actual <- matrix(m[rownames(ec),], ncol=1)
  expected.counts.df$md.var <- matrix(var, ncol=1)

  #add n
  background.dist.tmp <- background.dist[background.dist$peak %in% rownames(ec),]
  tmp <-  gene.sum[background.dist.tmp$gene,]
  expected.counts.df$n <- matrix(tmp, ncol=1) #this is breaking

  #get background distribution
  ec.background <- ec[,background.cells]
  ec.background <- matrix(ec.background, ncol=1)
  ec.background <- data.frame(expected.counts = ec.background)
  ec.background$peak <- rownames(ec)
  ec.background$actual <- matrix(m.background[rownames(ec),], ncol=1)
  ec.background$md.var <- matrix(var[,background.cells], ncol=1)

  tmp <-  gene.sum[background.dist.tmp$gene,background.cells]
  ec.background$n <- matrix(tmp, ncol=1)

  ec.background.sub <- ec.background[ec.background$n>0,]
  cutoff.ec <- quantile(ec.background.sub$expected.counts, 0.99)
  cutoff.n <- quantile(ec.background.sub$n, 0.99)
  max.n <- max(ec.background.sub$n[ec.background.sub$n < cutoff.n])
  max.ec <- max(ec.background.sub$expected.counts[ec.background.sub$expected.counts < cutoff.ec])

  min.n <- min(ec.background.sub$n[ec.background.sub$n < cutoff.n])
  min.ec <- min(ec.background.sub$expected.counts[ec.background.sub$expected.counts < cutoff.ec])

  lx <- number.bins
  ly <- number.bins
  n_step <- (max.n-min.n)/lx
  n_grid <- min.n + n_step*0:lx
  ec_step <- (max.ec-min.ec)/ly
  ec_grid <- min.ec + ec_step*0:ly

  tmp <- findInterval(ec.background.sub$expected.counts, ec_grid)
  tmp2 <-  findInterval(ec.background.sub$n, n_grid)

  ec.background.sub$ec.bin <- tmp
  ec.background.sub$n.bin <- tmp2

  ec.background.sub$ec_n <- paste0(ec.background.sub$ec.bin, "_", ec.background.sub$n.bin)
  df.sub <- ec.background.sub[unlist(lapply(split(1:nrow(ec.background.sub), ec.background.sub$ec_n),
                                            sample_within_groups, sample.n = sample.n)), ]
  x.matrix <- matrix(cbind(df.sub$n, df.sub$expected.counts), ncol=2)

  ### calculate regular grid for kernel estimates
  n_grid_midpoints <- calculate_midpoints(min.n, max.n, lx)
  ec_grid_midpoints <- calculate_midpoints(min.ec, max.ec, ly)

  grid <- matrix(cbind(rep(n_grid_midpoints, length(ec_grid_midpoints)),  rep(NA, length(n_grid_midpoints))), ncol=2)
  grid[,2] <- rep(ec_grid_midpoints, each=length(n_grid_midpoints))

  # calculate regularized estimates on regulatr grid
  mh <- kreg(x = x.matrix, y = df.sub$md.var, grid = grid)
  grid.var <- data.frame(reg.var = mh$y)
  grid.out <- mh$x
  colnames(grid.out) <- c("n", "ec")
  grid.var <- cbind(grid.var, grid.out)
  grid.var$n_bin <- rep(1:length(n_grid_midpoints), each = length(ec_grid_midpoints))
  grid.var$ec_bin <- rep(1:length(ec_grid_midpoints), length(n_grid_midpoints))
  grid.var$ec_n <- paste0(grid.var$ec_bin, "_", grid.var$n_bin)

  #now get variance in all data, not just NT
  expected.counts.df$ec.bin <- findInterval(expected.counts.df$expected.counts, ec_grid)
  expected.counts.df$n.bin <- findInterval(expected.counts.df$n, n_grid)
  expected.counts.df$ec_n <- paste0(expected.counts.df$ec.bin, "_", expected.counts.df$n.bin)

  t3 <- left_join(expected.counts.df, grid.var, by="ec_n")
  t3$reg.var[is.na(t3$reg.var)] <- t3$md.var[is.na(t3$reg.var)]
  t3$reg.var.new <- t3$reg.var
  t3$reg.var.new[t3$reg.var< min.variance] <- min.variance # variance threshold
  var.fit <- matrix(t3$reg.var.new, nrow=nrow(ec), ncol=ncol(ec))
  return(var.fit)
}


#generate midpoints
calculate_midpoints <- function(a, b, n_intervals) {
  width <- (b - a) / n_intervals
  first_point <- a + width/2
  midpoints <- seq(first_point, by = width, length.out = n_intervals)
  return(midpoints)
}


sample_within_groups <- function(x, sample.n) {
  if (length(x) <=  sample.n) return(x)
  x[x %in% sample(x, sample.n)]
}

