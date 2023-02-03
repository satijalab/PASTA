# Overview 

PASTA is a toolkit for analysis, visualization, and exploration of cellular heterogeneity in alternative polyadenylation (APA) from scRNA-seq datasets. While scRNA-seq has been widely applied to characterize cellular heterogeneity in gene expression levels, these datasets can also be leveraged to explore variation in transcript structure as well. For example, reads from 3’ scRNA-seq datasets can be leveraged to quantify the relative usage of multiple polyadenylation (polyA) sites within the same gene. 


In [Kowalski*, Wessels*, Linder*](link) et al, we introduce a statistical framework, based on the Dirichlet multinomial distribution, to characterize heterogeneity in the relative use of polyA sites at single-cell resolution. This framework enables upervised differential analysis (differential polyadenylation between groups of cells), but also unsupervised analysis and visualization.

We have implemented these methods in an R package PASTA (PolyA Site analysis using relative Transcript Abundance) that interfaces directly with Seurat. 	This repository represents an initial release of PASTA, we will be adding significant functionality and documentation in future releases.

# Installation

To install PASTA, please run:

```{R}
remotes::install_git("satijalab/PASTA")
```

# Vignette

Check out the [PASTA vignette](http://www.satijalab.org/seurat/pasta_vignette.html) to perform supervised and unsupervised analysis on an scRNA-seq dataset of human PBMC

# Links

[Seurat](www.satijalab.org/seurat): R package for the analysis, integration, and exploration of scRNA-seq data

[CPA-Perturb-seq](link): Preprint describing analytical methods for characterizing APA in single-cell datasets
