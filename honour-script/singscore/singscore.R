# singscore analysis 
# setup 
## for running on cluster
.Library <- "/home/ygao/R_LIB"
.libPaths(.Library)

if (!requireNamespace("singscore")) {
  library(BiocManager)
  BiocManager::install("singscore")
}

library(singscore)
source("../GSEA/parsePathways.R")

classPathDir <- function(ptw) {
  # ptw with _UP and _DN as suffix to indicate the directions 
  up <- grep("_UP$", names(ptw), value = T)
  dn <- grep("_DN$", names(ptw), value = T)
  undir <- setdiff(names(ptw), c(up, dn))
  up.pheno <- str_replace(up, "_UP$", "")
  dn.pheno <- str_replace(dn, "_DN$", "")
  dir.pheno <- intersect(up.pheno, dn.pheno)
  unpair.up <- setdiff(up, paste0(dir.pheno, "_UP"))
  unpair.dn <- setdiff(dn, paste0(dir.pheno, "_DN"))
  lst <- list(dir.pheno = dir.pheno, unpair.dir = c(unpair.up, unpair.dn), undir = undir)
  print(sapply(lst, "length"))
  return(lst)
}

runAndWriteSingscore <- function(gs, mode, ptw, ptw.name, labels, rd = rank.data, isfilter = F, filename = NULL) {
  # this function runs singscore and writes the results
  # mode is one of "paired.dir", "unpaired.dir", "undir"
  
  if (is.null(filename)) {
    filename <- paste(Sys.Date(), mode, ptw.name, isfilter, "results.txt", sep = "_")
    filename <- paste0("out/", filename)
  }
  dir <- mode != "undir"
  cn <- c("TotalScore", "TotalDispersion", "phenotype", "GAMuT_ID", 'p.val', "Type")
  if (mode == "paired.dir") {
    
    cn <- c("TotalScore", "TotalDispersion", "UpScore", "UpDispersion", 
            "DownScore", "DownDispersion", "phenotype", "GAMuT_ID", "p.val", "Type")
  } 
  ncores <- 18
  
  # write out the header of the table
  cat("file name:", filename, "\n")
  write.table(matrix(cn, 1, length(cn)), filename, col.names = F, row.names = F, 
              quote = F, sep = "\t", append = F)
  
  pb <- txtProgressBar(min = 0, max = length(gs), initial = 0)
  i <- 0
  for (ph in gs) {
    if (mode == "paired.dir") {
      up.genes <- GeneSet(geneIdType = SymbolIdentifier(), geneIds = ptw[[paste0(ph, "_UP")]])
      down.genes <- GeneSet(geneIdType = SymbolIdentifier(), geneIds = ptw[[paste0(ph, "_DN")]])
      scoredf <- simpleScore(rd, upSet = up.genes, downSet = down.genes, 
                             knownDirection = T)
      permuteResult <- generateNull(upSet = up.genes, downSet = down.genes, 
                                    rankData = rank.data, centerScore = T,
                                    knownDirection = T, B = 1000, ncores = ncores, 
                                    seed = 1, useBPPARAM = NULL)
      
    } else {
      up.genes <- GeneSet(geneIdType = SymbolIdentifier(), geneIds = ptw[[ph]])
      scoredf <- simpleScore(rd, upSet = up.genes, knownDirection = dir)
      permuteResult <- generateNull(upSet = up.genes, rankData = rank.data, 
                                    centerScore = T, knownDirection = dir, B = 1000, 
                                    ncores = ncores, seed = 1, useBPPARAM = NULL)
      
    }
    pvals <- getPvals(permuteResult = permuteResult, scoredf = scoredf)
    scoredf$phenotype <- ph
    scoredf$GAMuT_ID <- rownames(scoredf)
    scoredf$p.val <- pvals[scoredf$GAMuT_ID]
    scoredf$Type <- as.character(labels$Labels)
    write.table(scoredf, filename, col.names = F, row.names = F, 
                quote = F, sep = "\t", append = T)
    i <- i + 1
    setTxtProgressBar(pb, i)
    # for testing only
    # if (length(scorelist) > 5) break
  }
  close(pb)
}

# read data files
load("../GoldDataUseWhenEverPossible/mut.and.exp.matrix.RData")
filter <- F
t.exp.mat <- t(mocTumourExpression)
if (filter) {
  eom <- read.table("EOMexclude.csv", sep = ",", header = F, 
                    stringsAsFactors = F)$V1
  eom <- paste0("GAMuT_", eom)    # the samples to exclude from MOC
  t.exp.mat <- t.exp.mat[, setdiff(colnames(t.exp.mat), eom)]
}
b.exp.mat <- t(benExpression)
exp.mat <- cbind(t.exp.mat, b.exp.mat)

# form gene sets from downloaded pathways
ptw.file <- "Data/msigdb_v6.2_files_to_download_locally/msigdb_v6.2_GMTs/c2.all.v6.2.symbols.gmt"
ptw <- parsePathway(ptw.file, info.col = F)

# make sure each gene set contains unique genes (other wise error when creating gene sets)
tmp <- sapply(ptw, function(e) {length(which(duplicated(e))) > 0})
ptw[tmp] <- lapply(ptw[tmp], "unique")
# trim the ptw so that each gene set only contains the genes in RNAseq
ptw <- trimIncomplete(ptw, panel.genes = rownames(exp.mat))
cat("ptw length:", length(ptw), "\n")
# 3615

ptw.list <- classPathDir(ptw)

# prepare to run singscore
labels <- data.frame(Labels = c(rep("MOC", ncol(t.exp.mat)), rep("BEN", ncol(b.exp.mat))), 
                     row.names = colnames(exp.mat))
se <- SummarizedExperiment::SummarizedExperiment(assays = SimpleList(cpm = exp.mat), colData = labels)
rank.data <- rankGenes(se)
# scorelist <- list()   # not using scorelist because takes too much memory and too slow


## run singscore for each phenotype
## and write the score data frame along running 
runAndWriteSingscore(gs = ptw.list$dir.pheno, mode = "paired.dir", labels = labels, rd = rank.data, 
                     ptw.name = "c2.all", isfilter = filter, ptw = ptw)
runAndWriteSingscore(gs = ptw.list$unpair.dir, mode = "unpaired.dir", labels = labels, rd = rank.data,
                     ptw.name = "c2.all", isfilter = filter, ptw = ptw)
runAndWriteSingscore(gs = ptw.list$undir, mode = "undir", labels = labels, rd = rank.data, 
                     ptw.name = "c2.all", isfilter = filter, ptw = ptw)


