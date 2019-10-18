###### run DriverNet ---------
.Library <- "/home/ygao/R_LIB/"
.libPaths(.Library)

biocDepend <- c("DriverNet")
for (dep in biocDepend) {
  if (!requireNamespace(dep)) {
    
    BiocManager::install(dep, ask = FALSE)
  }
  library(package = dep, character.only = TRUE)
}

print("Success: loaded DriverNet")
# used for commandline job
# use.prev.data <- commandArgs(TRUE)
DESTIN_FOLDER <- "~/R/"
# if use symmetric influence graph

load("../GoldDataUseWhenEverPossible/influence.graph.all.RData")
load("../GoldDataUseWhenEverPossible/mut.and.exp.matrix.RData")

outFolder <- "./out/"

ca125 <- paste0("GAMuT_", read.csv("../CA125exclude.csv", header = F, stringsAsFactors = F)[,1])
# add the EOM case to ca125 to be removed 
ca125 <- c(ca125, "GAMuT_22826")

influenceGraph <- as.matrix(ggiANDppi)
# cut mutMatrix, mocTumourExpression, and influenceGraphs to genes they share
# use symmetrical matrix for this
genes.union <- unique(c(colnames(influenceGraph), rownames(influenceGraph),
                        colnames(mutMatrix), colnames(mocTumourExpression)))
mutMatrix <- mutMatrix[setdiff(rownames(mutMatrix), ca125), which(colnames(mutMatrix) %in% genes.union)]
influenceGraph <- influenceGraph[which(rownames(influenceGraph) %in% genes.union),
                                 which(colnames(influenceGraph) %in% genes.union)]
mocTumourExpression <- mocTumourExpression[setdiff(rownames(mocTumourExpression), ca125), intersect(colnames(mocTumourExpression), genes.union)]
if (any(rownames(mutMatrix) != rownames(mocTumourExpression))) {
  stop("Rownames of mutMatrix and mocTumourExpression differs, should check name matching")
}

# see input data
mocOutlierMatrix <- DriverNet::getPatientOutlierMatrix(mocTumourExpression)

inputDataList <- list(mutation.matrix = mutMatrix,
                      outlier.matrix = mocOutlierMatrix,
                      influence.graph = influenceGraph)
for (mat in names(inputDataList)) {
  print(mat)
  print(inputDataList[[mat]][1:10, 1:10])
  cat('\n\n')
}

# The main function to compute drivers
print("running computeDrivers")
driversList = computeDrivers(patMutMatrix = mutMatrix, 
                             patOutMatrix = mocOutlierMatrix,
                             influenceGraph = influenceGraph, 
                             outputFolder=outFolder, printToConsole=TRUE)


save(list = c("driversList"), file = paste0(DESTIN_FOLDER, Sys.Date(), "_drivernet_driverList.RData"))

# random permute the gene labels to compute p-values
print("running randomDriversResult")
randomDriversResult = computeRandomizedResult(patMutMatrix=mutMatrix,
                                              patOutMatrix=mocOutlierMatrix, influenceGraph=influenceGraph,
                                              geneNameList=factor(colnames(mutMatrix)), outputFolder=outFolder, printToConsole=TRUE,
                                              numberOfRandomTests=50, weight=FALSE, purturbGraph=FALSE, purturbData=TRUE)
save(list = c("driversList", "drivers", "randomDriversResult"), file = paste0(DESTIN_FOLDER, Sys.Date(), "_drivernet_randRes.RData"))
# Summarize the results
print("running resultSummary")
res = resultSummary(mainResult = driversList, 
                    randResult = randomDriversResult, 
                    patMutMatrix = mutMatrix,
                    influenceGraph = influenceGraph, 
                    outputFolder = outFolder, 
                    printToConsole=TRUE)
print("saveing session info")
save(list = c("driversList", "drivers", "randomDriversResult", "res"), file = paste0(DESTIN_FOLDER, Sys.Date(), "_drivernet_res.RData"))
# save the complete image
save.image(file = paste0("out/", Sys.Date(), "_drivernet_FullImage.RData"))