######## generate mutation matrix and expression matrix for DriverNet --------------
.libPaths("/home/ygao/R_LIB/")
.Library <- "/home/ygao/R_LIB/"   # change default lib.loc

dependencies <- c("devtools", 
                  "snowfall", 
                  "maxstat", 
                  "Matrix",
                  "tidyverse")

for (pack in dependencies) {
  if (!requireNamespace(pack)) {
    install.packages(pack, repos = "https://cran.ms.unimelb.edu.au/", quiet = TRUE)
  }
  library(package = pack, character.only = TRUE)
}

print("Success: loaded all required packages")

findMatchingRows <- function(ids, df, exclude.samples = NULL, targetCol = 1) {
  #find rows in df that has column targetCol matching all elemets in ids
  matchingRows <- c()
  noMatchIds <- c()
  multipleMatchIds <- c()
  matches <- c()
  for (id in ids) {
    if (id %in% exclude.samples) {next}
    
    currentMatch <- df[grep(id, df[, targetCol], fixed = TRUE),]
    matchedId <- unique(as.vector(currentMatch[, targetCol]))
    #check for number of matches
    if (length(matchedId) == 0) {
      noMatchIds <- c(noMatchIds, id)
    }
    #check for multiple matches
    else if (length(matchedId) > 1) {
      multipleMatchIds <- c(multipleMatchIds, id)
      cat(id, "matches multiple ids:", matchedId, "\n")
      
    }
    #check if the match occured before
    else if (matchedId %in% matches) {
      cat(unique(as.vector(currentMatch[,targetCol])) , "\n")
    }
    else {
      cat(id, "matches", matchedId, "\n")
      matches <- c(matches, unique(as.vector(currentMatch[,targetCol])))
      currentMatch[,targetCol] <- rep(paste0("GAMuT_", id), nrow(currentMatch))
      matchingRows <- rbind(matchingRows, currentMatch)
    }
  }
  if (length(noMatchIds) > 0) {cat(noMatchIds, "found no match\n")}
  return(list(matchingRows=as.data.frame(matchingRows), noMatches=noMatchIds))
}

# set if CNV data should be filtered
highCNV <- TRUE

# load data from file
filenames <- list(expression = "mutExpData/expressionData.csv", 
                  variant = "mutExpData/mocVariantData.csv", 
                  cnv = "mutExpData/mocCNVData.txt")

# set ids of MOC samples
mocGamutIDs <- c('GAMuT_1000', 'GAMuT_1818', 'GAMuT_20077', 'GAMuT_20104', 'GAMuT_22826', 'GAMuT_23022', 'GAMuT_23836', 'GAMuT_32119', 'GAMuT_34951', 'GAMuT_51107', 'GAMuT_51149', 'GAMuT_51163', 'GAMuT_60259', 'GAMuT_134016', 'GAMuT_134032', 'GAMuT_134035', 'GAMuT_51107R', 'GAMuT_844.07G3', 'GAMuT_HOV449', 'GAMuT_IC138', 'GAMuT_IC257', 'GAMuT_OV0711', 'GAMuT_OV0772', 'GAMuT_OV0841', 'GAMuT_OV1217', 'GAMuT_OV1417', 'GAMuT_OV1621', 'GAMuT_OV1723', 'GAMuT_OV1754', 'GAMuT_OV1754R', 'GAMuT_OV1991', 'GAMuT_OV2090', 'GAMuT_OV2405', 'GAMuT_OV2602.A', 'GAMuT_UNSW135425', 'GAMuT_VOA1675', 'GAMuT_WM1070A', 'GAMuT_WM1250A', 'GAMuT_WM1339A', 'GAMuT_WM1445A', 'GAMuT_WM1747A', 'GAMuT_WM982A', 'GAMuT_WM982B')
#mocIDs without the "GAMuT_" prefix for mutation data files
mocIDs <- c('1000', '1818', '20077', '20104', '22826', '23022', '23836', '32119', '34951', '51107', '51149', '51163', '60259', '134016', '134032', '134035', '51107R', '844.07G3', 'HOV449', 'IC138', 'IC257', 'OV0711', 'OV0772', 'OV0841', 'OV1217', 'OV1417', 'OV1621', 'OV1723', 'OV1754', 'OV1754R', 'OV1991', 'OV2090', 'OV2405', 'OV2602.A', 'UNSW135425', 'VOA1675', 'WM1070', 'WM1250A', 'WM1339A', 'WM1445A', 'WM1747A', 'WM982A', 'WM982B')
noData <- c("134032", "51107R", "OV1754R", "UNSW135425", "WM1250A", "WM1339A", "WM1445A", "WM1747A", "WM982A")

# variant data -------
# read in variant file
variantData <- read.csv(filenames$variant, as.is = T)
variantData <- variantData[grep("^[[:digit:]]*$", variantData$Consequence_Rank), ]
variantData$Consequence_Rank <- as.numeric(variantData$Consequence_Rank)
variantData <- variantData[variantData$Consequence_Rank < 5, ]

#search for MOC cases for variant data
resultVariant <- findMatchingRows(ids = mocIDs, df = variantData, exclude.samples = noData, targetCol = 1)
mocVariantData <- resultVariant$matchingRows
mocVariantNoMatches <- resultVariant$noMatches

# CNV Data -------
CNVData <- read.delim(filenames$cnv, as.is = T)
if (highCNV) CNVData <- CNVData[CNVData$Event %in% c("Homozygous Copy Loss", "High Copy Gain"),]
resultCNV <- findMatchingRows(mocIDs, CNVData, exclude.samples = noData, targetCol = 1)
mocCNVData <- resultCNV$matchingRows
mocCNVNoMatches <- resultCNV$noMatches

#construct mutation matrix ------
#unionSampleIds are the row names, colnames are the gene names
# Unsure whether or not to use UNION,  maybe Intersection?
unionSampleIds <- union(mocCNVData$Sample, mocVariantData$Sample)
mutMatrix <- matrix(nrow = length(unionSampleIds), ncol = 0)
rownames(mutMatrix) <- unionSampleIds

# running extremely slow
# try sparse matrix
# try collect gene names -> construct matrix with fixed size
# -> coerce sparse matrix -> iterate thru variant and CNV data
# -> change numbers

# get gene names
genesMutated <- unique(c(as.vector(mocVariantData$SYMBOL), unlist(strsplit(as.vector(mocCNVData$Gene.Symbols), ", "))))

#construct matrix with "known" size
mutMatrix <- matrix(0, nrow = length(unionSampleIds), 
                    ncol = length(genesMutated), 
                    dimnames = list(unionSampleIds, genesMutated))
mutMatrix <- as(mutMatrix, "sparseMatrix")

for (id in unionSampleIds) {
  
  #collect Variant data
  if (id %in% mocVariantData$Sample) {
    mutMatrix[id, as.vector(mocVariantData$SYMBOL[mocVariantData$Sample == id])] = 1
  }
  #get CNV data
  if (id %in% mocCNVData$Sample) {
    mutMatrix[id, unlist(strsplit(as.vector(mocCNVData$Gene.Symbols[mocCNVData$Sample == id]), ", "))] = 1
  }
}
#sort columns by gene names
mutMatrix <- as.matrix(mutMatrix[, order(colnames(mutMatrix))])
cat(length(unionSampleIds), "samples with", length(genesMutated), "genes included in mutMatrix\n")
# mutMatrix is the binary patient mutation matrix
#comment next line to disable exporting mutMatrix
write.csv(mutMatrix, file = "mutationMatrixPatientByGene.csv", quote = FALSE, row.names = TRUE)
print("Success: Patient mutation matrix built")
# expression data

tumourExpression <- read.csv(filenames$expression, header = TRUE)

# check duplicated gene names, and check if rows are exact duplications
dupGenesTumour <- unique(tumourExpression$X[duplicated(tumourExpression$X)])
cat(length(dupGenesTumour), "duplicated genes found:", dupGenesTumour, "\n")

# assign gene names to rownames and delete the gene name column
rownames(tumourExpression) <- tumourExpression$X
tumourExpression$X <- NULL

# delete experimental repeated rows ()
repeats <- c("GAMuT_51107", "GAMuT_OV1754")
# GAMuT_51107 t.test significant - original mean: 4.362801 repeat: 4.578879 
# GAMuT_OV1754 t.test non significant = original mean: 4.430244 repeat: 4.452958
# will delete orignal rows for now (assume repeats more accurate)
tumourExpression[, repeats] <- NULL
tumourExpression <- as.data.frame(t(tumourExpression))

print("Success: loaded tumourExpressionMatrix")

# get MOC samples only
tumourExpression = cbind(rownames(tumourExpression), tumourExpression)
resultExpression <- findMatchingRows(unionSampleIds, tumourExpression, exclude.samples = noData)
mocTumourExpression <- resultExpression$matchingRows[,-1]
expressionNoMatches <- resultExpression$noMatches

# change the rownames of mocTumourExpression because the names are different
rn <- rownames(mocTumourExpression)
rn[rn == "GAMuT_51107R"] <- "GAMuT_51107"
rn[rn == "GAMuT_OV1754R"] <- "GAMuT_OV1754"
rn[rn == "GAMuT_WM1070A"] <- "GAMuT_WM1070"
rownames(mocTumourExpression) <- rn

print("Success: generated mutation and expression matrix")

# DESTIN_FOLDER <- "~/R/"

# to save
# save.image(paste0(DESTIN_FOLDER, Sys.Date(), "getInput.RData"))

# save(list = c("mutMatrix",
#               "mocCNVData",
#               "mocVariantData",
#               "mocTumourExpression",
#               "mocCNVNoMatches",
#               "expressionNoMatches",
#               "mocVariantNoMatches",
#               "highCNV"), file = "out/mut.and.exp.matrix.RData")
# 

# save(list = c("mutMatrix",
#               "mocCNVData",
#               "mocVariantData",
#               "mocTumourExpression",
#               "mocCNVNoMatches",
#               "expressionNoMatches",
#               "mocVariantNoMatches",
#               "highCNV"), file = "../GoldDataUseWhenEverPossible/mut.and.exp.matrix.RData")
# print("Success: data saved")
