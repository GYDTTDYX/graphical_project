##
## CStaR Script File
## =================
##
## Description:
## ------------
## Script for running CStaR on a data set of Arabidopsis thaliana as
## as described in "Causal Stability Ranking", Stekhoven et al., 2012,
## not published.
##
## Author: Daniel J. Stekhoven, stekhoven@stat.math.ethz.ch
###############################################################################


## Dependencies on external packages ==========================================
library(pcalg) # IDA
## ============================================================================


## Input parameters ===========================================================
alpha <- 0.1                 # tuning parameter for PC algorithm
noStabRuns <- 100              # number of stability runs (default=100)
corCutoff <- 0.05             # cutoff for marginal correlation w.r.t. response
exprCutoff <- 4.5            # log2 gene expression cutoff
qSet <- seq(2000, 100, -50)  # set of 'q' values used for summary ranking
## ============================================================================


## DATA =======================================================================
dat <- read.table('final_count_data_100_stage_drop_na.tsv', header=TRUE, sep="\t")
dat <- dat[c("Outcome", setdiff(names(dat), "Outcome"))]
## Pre-processing -------------------------------------------------------------

yInd <- 1                           # index of response
## Gene names -----------------------------------------------------------------
## ============================================================================

## CStaR PROCEDURE ============================================================
t.start <- proc.time()                         # start timer
t.date <- format(Sys.time(), "%Y_%m_%d_%H.%M") # time stamp for file saving

## INITIATE STABILITY SELECTION -----------------------------------------------
n <- nrow(dat)
subSize <- floor(n/2)                 # size of subsample
p <- ncol(dat)-1                      # number of covariables
effList <- vector('list', noStabRuns) # list of resulting stabilities
effMat <- matrix(0, ncol=p, nrow=noStabRuns) # table of causal effects
colnames(effMat) <- colnames(dat)[-1]

## STABILITY SELECTION --------------------------------------------------------
for (k in 1:noStabRuns){
  subIndices <- sample(n, subSize, replace = FALSE)
  subDat <- dat[subIndices,] # subsample the data set
  ## Pre-select data w.r.t. marginal correlation with the response
  marCor <- apply(subDat, 2, function(x) abs(cor(subDat[,yInd], x)))
  subDat <- subDat[,marCor>=corCutoff]
  subP <- ncol(subDat)       # Note: p+1 (the response)
  subDat <- scale(subDat)    # standardize to mean zero and sd one
  effects <- numeric(subP-1) # vector of causal effects

  ## IDA ----------------------------------------------------------------------
  ## estimate equivalence class (CPDAG) on subsample
  indepTest <- gaussCItest
  suffStat <- list(C=cor(subDat), n=subSize)
  pcFit <- pc(suffStat, indepTest, p=subP, alpha, verbose=FALSE)
  
  ## estimate lower bound for causal effects ----------------------------------
  for (j in 2:(subP)){
    multSet <- ida(j, yInd, cov(subDat), pcFit@graph, method='local')
    effMat[k,colnames(subDat)[j]] <- min(abs(multSet))  
  }
} # end of for (k in 1:noStabRuns)

## SUMMARY RANKING ------------------------------------------------------------
medianEffect <- apply(effMat, 2, median) # median causal effect
maxExpr <- apply(dat[,-1], 2, max) # maximum log2 gene expression
P <- length(qSet)
dimnames_list <- list(
  stabRuns = paste0("run", 1:noStabRuns),  # Names for the first dimension
  genes = paste0("gene", 1:p),            # Names for the second dimension
  samples = paste0("sample", 1:P)         # Names for the third dimension
)
topArr <- array(FALSE, c(noStabRuns, p, P), dimnames_list)
stabMat <- matrix(0, ncol=p, nrow=P) # matrix to store frequencies
colnames(stabMat) <- colnames(dat[,-1])
rownames(stabMat) <- qSet
rankMat <- stabMat # matrix to store ranks
errorMat <- stabMat # matrix to store PCER

PCER <- function(p, q, freq=NA){
  ## Purpose:
  ## Calculate per-comparison error rate (PCER) given a selection frequency.
  ## ----------------------------------------------------------------------
  if (all(is.na(freq))){
    return(NULL)
  } else {
    PCER <- q^2 / (2*freq - 1) / p^2
  }  
}
## Calculate for each gene its frequency, rank and error for different q's
for (t.q in 1:P){
  for (t.run in 1:noStabRuns){
    lowerQbound <- sort(effMat[t.run,], decreasing=TRUE)[t.q]
    if (lowerQbound == 0){
      lowerQbound <- min(effMat[t.run, effMat[t.run,] != 0])
    }
    t.include <- effMat[t.run,] >= lowerQbound
    topArr[t.run,,t.q] <- t.include
  }
  stabMat[t.q,] <- apply(topArr[,,t.q], 2, function(x) sum(x)/n)
  rankMat[t.q,] <- rank(-stabMat[t.q,], ties.method = "min")
  errorMat[t.q,] <- PCER(p=p, q=qSet[t.q], freq=stabMat[t.q,])
}
errorMat[errorMat<0] <- 1 # replace negative errors with PCER == 1
# --------------------------------------------------------------------------
medianRank <- apply(rankMat, 2, median); names(medianRank) <- colnames(rankMat)
medianPCER <- apply(errorMat, 2, median)
# --------------------------------------------------------------------------
stabRank <- data.frame(medianRank, medianEffect, maxExpr, medianPCER)#, name=annotation[,3]) # generate causal stable ranking
stabRank <- stabRank[order(stabRank[,1], -stabRank[,2]),]

library(biomaRt)
stabRank$Gene = rownames(stabRank)
genes = stabRank$Gene
genes_no_version <- sub("\\..*", "", genes)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",                     
  values = genes_no_version,                       
  mart = ensembl
)
stabRank$Gene_Symbol <- gene_mapping$hgnc_symbol[match(sub("\\..*", "", stabRank$Gene), gene_mapping$ensembl_gene_id)]
write.csv(stabRank, "gene_analysis_results_Cstar_metastasis.csv", row.names = FALSE)
