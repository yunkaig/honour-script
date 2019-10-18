###### run DawnRank ---------

.Library <- "/home/ygao/R_LIB/"   # change default lib.loc
.libPaths("/home/ygao/R_LIB/")

library(devtools)
if (!requireNamespace("DawnRank")) {
  install_github("MartinFXP/DawnRank", quiet = TRUE)
}
library(DawnRank)
library(maxstat)
library(readr)

print("Success: loaded DawnRank")
# used for commandline job
# use.prev.data <- commandArgs(TRUE)
DESTIN_FOLDER <- "out/"
# if use symmetric influence graph
use.symmetric <- TRUE

load("../GoldDataUseWhenEverPossible/influence.graph.all.RData")
load("../GoldDataUseWhenEverPossible/all.variant.and.cnv.RData")
load("../GoldDataUseWhenEverPossible/mut.and.exp.matrix.RData")

influenceGraph <- as.matrix(ggiANDppi)


EOM <- c("GAMuT_22826")

# cut mutMatrix, mocTumourExpression and influenceGraphs to genes they share

# use symmetrical matrix for this
genes.all.inf.graph <- intersect(union(colnames(influenceGraph), rownames(influenceGraph)),
                                 intersect(colnames(mutMatrix), colnames(mocTumourExpression)))
mutMatrix <- mutMatrix[setdiff(rownames(mutMatrix), EOM), genes.all.inf.graph]
mocTumourExpression <- mocTumourExpression[setdiff(rownames(mocTumourExpression), EOM), genes.all.inf.graph]
benExpression <- benExpression[, genes.all.inf.graph]
influenceGraph <- influenceGraph[genes.all.inf.graph,
                                 genes.all.inf.graph]

# transpose the matrices for DawnRank
mutMatrix <- t(mutMatrix)
mocTumourExpression <- t(mocTumourExpression)
benExpression <- t(benExpression)

# see input data

inputDataList <- list(mutation.matrix = mutMatrix,
                      tumour.expresion.matrix = mocTumourExpression,
                      ben.expression.matrix = benExpression,
                      influence.graph = influenceGraph)
for (mat in names(inputDataList)) {
  print(mat)
  print(inputDataList[[mat]][1:10, 1:10])
  cat('\n\n')
}
# DawnRank main algorithm ------
# set gold standards (known driver genes)
knownDrivers <- c("KRAS", "TP53")

# normalize the tumor and normal data to get the differential expression
normalizedDawn <- DawnNormalize(tumorMat = mocTumourExpression, 
                                normalMat = benExpression)

# get the DawnRank Score Get some coffee, this might take a while!
dawnRankScore <- DawnRank(adjMatrix = influenceGraph, 
                          mutationMatrix = mutMatrix, 
                          expressionMatrix = normalizedDawn, 
                          mu = 3, 
                          goldStandard = knownDrivers, 
                          parallel = 8)

# look at the DawnRank scores for a few patients
dawnRankFrame <- dawnRankScore[[3]]
head(dawnRankFrame)

# get the aggregate DawnRank scores Get some coffee, this might take a
# while!
aggregateDawnRankScore <- condorcetRanking(scoreMatrix = dawnRankScore[[2]], 
                                           mutationMatrix = mutMatrix, 
                                           parallel = 8)

# look at top 10 ranked genes
top10 <- aggregateDawnRankScore[[2]][1:10]
top10

# get the individual cutoff for patient
dawnRankFrame$isCGC <- dawnRankFrame$isGoldStandard
out.dir <- paste0(Sys.Date(), "_dawnRankOut")
if (!(out.dir %in% dir())) {
  dir.create(out.dir)
}

#NOTE: the latest version of mvnorm should be installed
ranked.genes <- list()
all.results <- NULL
############13 pat is wrong, look at this 
###############
for (pat in colnames(mocTumourExpression)) {
  significance <- patspeccutoff(patient = pat, ms = dawnRankFrame, 
                              default = 95)
  # look for signficance. 
  # significance[[1]] contains the rows of dawnRankFrame for an individual patient
  # plus a binary column named "significant"
  
  # write tab-delim files for each patient
  write_delim(significance[[1]], paste0(out.dir, "/", pat, "_significance.txt"), delim = '\t')
  
  # pick out the ranked genes only and write out to one single file
  # to allow patient-wise comparisons
  ranked.genes[[pat]] <- as.character(significance[[1]][significance[[1]]$significant == 1, "Gene"])
  all.results <- rbind(all.results, significance[[1]])
}
max.len <- max(sapply(ranked.genes, length))
ranked.genes <- lapply(ranked.genes, FUN = function(e) {
  return(c(e, rep("", max.len - length(e))))
})
ranked.genes.df <- as.data.frame(ranked.genes, stringsAsFactors = FALSE)
write_delim(ranked.genes.df, paste0(out.dir, "/ranked.genes.all.txt"), delim = '\t')
write_delim(all.results, paste0(out.dir, "/all.results.txt"), delim = "\t")

# save the complete image
filename <- paste0(out.dir, "_DawnRank_FullImage.RData")
save.image(file = filename)
cat("Image saved to", filename, "\n")

#### problem pat: GAMuT_HOV449, GAMuT_WM982B
# warnings in running patspeccutoff
# In pmaxstat(STATISTIC, y, m) : value out of range in 'gammafn'
