
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
#' @param do.center Return the centered residuals. Default is TRUE.
#' @param do.scale Return the scaled residuals. Default is TRUE.
#' @param residuals.max Clip residuals above this value. Default is NULL (no clipping).
#' @param residuals.min Clip residuals below this value. Default is NULL (no clipping).
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
                               gene.names = "symbol",
                               min.counts.background = 5,
                               min.variance = 0.1,
                               do.scale = FALSE,
                               do.center = FALSE,
                               residuals.max = NULL, 
                               residuals.min = NULL,
                               verbose=TRUE)
 {
  if(verbose) {
    message("Calculating background distribution")
  }

  #if features in NULL, then specify all features in polyA assay
  #TO DO: change to all features with a gene annotation
  if (is.null(features)) {
    features <- rownames(GetAssayData(object, assay=assay))
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
    message(paste0("Using", background, " as background distribution"))
  }


  #check if symbols are contained in meta features
  if (!(gene.names %in% colnames(object[[assay]][[]]))) {
    stop("Gene.names column not found in meta.features, please make sure 
         you are specific gene.names correctly")
  }
  
  if (sum(is.na(object[[assay]][[]][features,gene.names])) > 0) {
    features.no.anno <- features[is.na(object[[assay]][[]][features,gene.names])]
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
  features.use <- background.dist$peak

  ##############################################################################
  #calculate sum of counts within each gene
  m <- GetAssayData(object = object, assay=assay, slot="counts")
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

  ec <- res[[1]]$ec
  var <- res[[1]]$var
  for(i in 2:length(res)) {
    if (!(is.null(res[[i]]))) {
      ec <- cbind(ec, res[[i]]$ec)
      var <- cbind(var, res[[i]]$var)
    }
  }

  ec <- t(ec)
  var <- t(var)
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
                      min.variance = min.variance)

  #calculate residual matrix
  residual.matrix <- (m-ec) / sqrt(var.reg)
  residual.matrix <- as.matrix(residual.matrix, nrow = nrow(residual.matrix))
  #M1 <- as(residual.matrix, "dgCMatrix")
  
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
  object[[assay]]@scale.data <- residual.matrix
  object <- LogSeuratCommand(object = object)
  return(object)
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
  nt.pseudo$gene <- paste0(object[[assay]][[]][features, gene.names], "_", object[[assay]]@meta.features[features, "strand"])
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
#'
#' @return Returns a list where first element is matrix of expected values for each peak within the genes,
#' second value is matrix of variance for each peak within the gene
#'
#' @importFrom MGLM MGLMfit
#' @concept residuals
#'
DirichletMultionmial <- function(gene.test, background.dist, m.background, gene.sum, ncells = ncells) {
  peaks <- background.dist$peak[background.dist$gene == gene.test]
  t <- m.background[rownames(m.background) %in% peaks,]
  t <- t(t)
  fit <- try(compareFit <- suppressWarnings(MGLMfit(t, dist="DM")), silent=TRUE)

  if (class(fit)!="try-error") {
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
#' @concept residuals
#'
#'
RegDMVar <- function(ec, var, m, m.background, background.dist, gene.sum, background.cells, min.variance = min.variance) {
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

  lx <- ly <- 30
  n_step <- (max.n-min.n)/lx
  n_grid <- min.n + n_step*0:lx
  ec_step <- (max.ec-min.ec)/ly
  ec_grid <- min.ec + ec_step*0:ly

  tmp <- findInterval(ec.background.sub$expected.counts, ec_grid)
  tmp2 <-  findInterval(ec.background.sub$n, n_grid)

  ec.background.sub$ec.bin <- tmp
  ec.background.sub$n.bin <- tmp2

  sample_n  = 1000

  #df.sub <- ec.background.sub %>% group_by(ec.bin, n.bin) %>% slice_sample(n = sample_n) #dplyr solution

  sample_within_groups <- function(x, sample_n) {
    if (length(x) <=  sample_n) return(x)
    x[x %in% sample(x, sample_n)]
  }
  ec.background.sub$ec_n <- paste0(ec.background.sub$ec.bin, "_", ec.background.sub$n.bin)
  df.sub <- ec.background.sub[unlist(lapply(split(1:nrow(ec.background.sub), ec.background.sub$ec_n), sample_within_groups, sample_n = sample_n)), ]
  x.matrix <- matrix(cbind(df.sub$n, df.sub$expected.counts), ncol=2)

  #kernel regularization
  pcf <- pcf.kernesti(x.matrix, y= df.sub$md.var ,h=0.5, N=c(lx, ly), kernel="gauss")
  df <- data.frame(cbind(pcf$value, pcf$down, pcf$high))
  colnames(df) <- c("reg.var", "n_low", "ec_low", "n.bin", "ec.bin")
  df$ec_n <- paste0(df$n.bin, "_", df$ec.bin)

  #now get variance in all data, not just NT
  expected.counts.df$ec.bin <- findInterval(expected.counts.df$expected.counts, ec_grid)
  expected.counts.df$n.bin <- findInterval(expected.counts.df$n, n_grid)

  #this is slow if you don't use dplyr left merge
  t3 <- left_join(expected.counts.df, df, by=c("n.bin", "ec.bin"))

  #if outside of 99th percentile, just use estimated variance
  #set minimum threshold to 0.1
  t3$reg.var[is.na(t3$reg.var)] <- t3$md.var[is.na(t3$reg.var)]
  t3$reg.var.new <- t3$reg.var
  t3$reg.var.new[t3$reg.var< min.variance] <- min.variance #set threshold as 0.1

  #filter(t3, actual==0 & expected.counts.nt>0) %>% head()
  var.fit <- matrix(t3$reg.var.new, nrow=nrow(ec), ncol=ncol(ec))

  return(var.fit)
}


#should we re-write this becuase this package is no longer on CRAN?
pcf.kernesti<-function(x,y,h,N,kernel="gauss",support=NULL)
  #from here https://github.com/cran/regpro/blob/master/R/pcf.kernesti.R#L77
{
  d = 2
  if (kernel=="bart")
    ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
  if (kernel=="gauss")
    ker<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
  if (kernel=="uniform")
    ker<-function(xx){ return( (rowSums(xx^2) <= 1) ) }

  if (kernel=="gauss") radi<-2*h else radi<-h

  recnum<-prod(N)
  value<-matrix(0,recnum,1)
  index<-matrix(0,recnum,2)

  support<-matrix(0,4,1)
  for (i in 1:2){
    support[2*i-1]<-min(x[,i])
    support[2*i]<-max(x[,i])
  }
  lowsuppo<-matrix(0,2,1)
  for (i in 1:2) lowsuppo[i]<-support[2*i-1]
  step<-matrix(0,2,1)
  for (i in 1:2) step[i]<-(support[2*i]-support[2*i-1])/N[i]

  for (i in 1:recnum){
    inde<-digit(i-1,N)+1
    arg<-lowsuppo+step*inde-step/2
    argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
    neigh<-(rowSums((argu-x)^2) <= radi^2)
    if (sum(neigh)>=2){     # if there are obs in the neigborhood

      xred<-x[neigh,]
      yred<-y[neigh]
      argu<-matrix(arg,dim(xred)[1],d,byrow=TRUE)

      w<-ker((xred-argu)/h)/h^d
      w<-w/sum(w)
      valli<-w%*%yred
    }
    else valli<-mean(y)

    value[i]<-valli
    index[i,]<-inde
  }

  down<-index-1
  high<-index

  pcf<-list(
    value=value,index=index,
    down=down,high=high,
    support=support,N=N)

  return(pcf)
}


digit<-function(luku,base){
  #Gives representation of luku for system with base
  #
  #luku is a natural number >=0
  #base is d-vector of integers >=2, d>=2,
  #base[d] tarvitaan vain tarkistamaan onko luku rajoissa
  #
  #Returns d-vector of integers.
  #
  #example: digit(52,c(10,10)), returns vector (2,5)
  #
  d<-length(base)
  digi<-matrix(0,d,1)
  jako<-matrix(0,d,1)
  jako[d]<-base[1]
  for (i in (d-1):1){
    jako[i]<-base[d-i+1]*jako[i+1]
  }
  vah<-0
  for (i in 1:(d-1)){
    digi[i]<-floor((luku-vah)/jako[i+1]) #if digi[i]>base[i], then ERROR
    vah<-vah+digi[i]*jako[i+1]
  }
  digi[d]<-luku-vah
  # annetaan vastaus kaanteisesti se 2354 annetaan c(4,5,3,2)
  # talloin vastaavuus sailyy base:n kanssa
  #apu<-matrix(0,d,1)
  #for (i in 1:d){
  #  apu[i]<-digi[d-i+1]
  #}
  apu<-digi[d:1]
  return(apu)
}






