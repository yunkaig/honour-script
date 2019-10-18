# create personal library and load packages

dependencies <- c("BiocManager", 
                  "devtools", 
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

########### START OF "generateInfluenceGraph.R" ###################
# new method: create influence graph from local datasets ----
# helper functions -----
to.char.data.frame <- function(df) {
  for (i in 1:length(df[1, ])) df[, i] <- as.character(df[, i])
  return(df)
}

geneSymbolFromID <- function(convDF, ids) {
  # convDF: first row - gene ids; second row - gene symbols
  ids <- as.vector(ids, mode = "character")
  convDF <- distinct(convDF)
  convDF[,1] <- as.vector(convDF[,1], mode = "character")
  convDF[,2] <- as.vector(convDF[,2], mode = "character")
  
  gene.symbols <- vapply(ids, FUN = function(x) {
    if (is.na(x)) return(as.character(NA))
    if (x %in% convDF[,1]) {
      return(paste(convDF[which(convDF[,1] == x), 2], collapse = "|"))
    }
    return(as.character(NA))
  }, FUN.VALUE = character(1), USE.NAMES = FALSE)
  return(gene.symbols)
}

symbolForDatabase <- function(ids, use.internal = TRUE) {
  # this function takes in a data.frame of ids, single type of id 
  # within each column (the output of geneSymbolFromID),
  # and returns a data.frame with the same dimensions
  # but converted to symbols
  
  sym.data.frame <- data.frame(row.names = 1:length(ids[,1]))
  if (use.internal) convDF <- mappingFileList$hgnc.map.internal
  else convDF <- mappingFileList$hgnc.map.external
  for (col.name in colnames(ids)) {
    curr.column <- ids[, col.name]
    id.type <- strsplit(col.name, "\\.")[[1]][1]
    curr.convDF <- convDF[, mappingCols$hgnc.map[c(id.type, "symbol")]]
    sym.data.frame <- cbind(sym.data.frame, data.frame(X = geneSymbolFromID(curr.convDF, curr.column), stringsAsFactors = FALSE))
    
  }
  colnames(sym.data.frame) <- paste0("symbol.from.", colnames(ids))
  sym.data.frame <- cbind(sym.data.frame, symbol = apply(sym.data.frame, MARGIN = 1, FUN = function(x) return(paste(unique(as.character(x[!is.na(x)])), collapse = "|"))))
  return(sym.data.frame)
}

getIDtypes <- function(ids, keys) {
  # ids is a data frame of an arbituary number of columns 
  # this function will return the likely names of idtypes in this data set
  idTypeKey <- list()
  for (i in 1:length(ids[, 1])) {
    break.down <- unlist(strsplit(as.character(unlist(ids[i, ], use.names = FALSE)), ":|\\|"), use.names = FALSE)
    for (key in keys) {
      selected.ids <- unique(grep(key, break.down, ignore.case = TRUE, value = TRUE))
      if (!is_empty(selected.ids)) {
        if (length(selected.ids) > 1) warning(paste("Multiple matches were found for", key))
        idTypeKey[[key]] <- unique(c(selected.ids, idTypeKey[[key]]))
      }
    }
    
    if (all(keys %in% names(idTypeKey))) return(idTypeKey)
  }
  
  return(idTypeKey)
}

separateByIdtype <- function(id.columns, idTypeNames) {
  # this function takes in a data.frame of ids, and return
  # a dataframe with separated ids (id only, without prefix)
  sepDataFrame <- apply(id.columns, MARGIN = 1, FUN = function(x) {
    x <- unlist(strsplit(unlist(x, use.names = FALSE), "\\|"))
    sep.id.vector <- character(0)
    for (id.name in idTypeNames) {
      sep.ids <- sapply(strsplit(grep(id.name, x,value = TRUE), ":"), FUN = function(e) return(e[2]))
      if (length(sep.ids) == 0) sep.id.vector <- c(sep.id.vector, NA)
      else sep.id.vector <- c(sep.id.vector, paste(sep.ids, collapse = "|"))
    }
    return(sep.id.vector)
  })
  
  sepDataFrame <- data.frame(t(sepDataFrame), row.names = NULL, stringsAsFactors = FALSE)
  colnames(sepDataFrame) <- names(idTypeNames)
  return(sepDataFrame)
}

processInteractionData <- function(db) {
  # this function takes in the database name and process the data.frame 
  # if this database has symbols, SKIP steps 1-2, id.type = NA
  # 1. get what id type names were used (returned in a list)
  # 2. separate columns by id names 
  # 3. delete columns that were not used
  # 4. rename the columns: detection.method, evidence.type, interaction.type, source.database, and ids
  # returned value is a list with id type names and formated data.frame
  
  int.data <- interactionDataList[[db]]
  # data.cols: columns with accessory data
  # id.cols: columns with id data
  # not.used.cols : columns that will be deleted
  data.cols <- interactionDataCols[[db]][!(names(interactionDataCols[[db]]) %in% c("A", "B"))]
  id.cols <- interactionDataCols[[db]][c("A", "B")]
  not.used.cols <- c(1:length(int.data[1,]))[-unlist(interactionDataCols[[db]], use.names = FALSE)]
  
  if (db %in% names(db.with.symbols)) id.type.names <- NA 
  else {
    
    # find id type names used in the database
    id.type.names <- getIDtypes(int.data[, id.cols$A], keys = used.idtypes)
    id.types <- unlist(id.type.names)[used.idtypes[used.idtypes %in% names(id.type.names)]]
    # separate id columns by id types
    for (intr in c("A", "B")) {
      sep.data.frame <- separateByIdtype(int.data[, id.cols[[intr]]], id.type.names)
      colnames(sep.data.frame) <- paste(colnames(sep.data.frame), intr, sep = ".")
      int.data <- cbind(int.data, sep.data.frame)
    }
  }
  # delete rows
  colnames(int.data)[unlist(data.cols)] <- names(data.cols)
  int.data[, not.used.cols] <- NULL
  for (i in 1:length(int.data[1,])) int.data[,i] <- as.character(int.data[,i])
  return(list(id.type = id.type.names, new.data.frame = int.data))
}

# to save data.list
save.new.int.data <- function(data.list, newIntDataFolder = processed.data.folder) {
  data.name <- deparse(substitute(data.list))
  for (db in names(data.list)) {
    file.name <- paste0(newIntDataFolder, Sys.Date(), "_", db, "_", data.name, ".txt")
    print(file.name)
    write.table(data.list[[db]], file.name, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
}

interactionDataFolder <- "interactionData/"
# load(".RData")
# loadhistory(".Rhistory")

# sources <- read_file(paste0(interactionDataFolder, "sources.txt"))

interactionFileLocations <- list(reactome = paste0(interactionDataFolder, "reactome.interactions.txt"),
                                 innatedb = paste0(interactionDataFolder, "innate.interactions.txt"),
                                 biogrid = paste0(interactionDataFolder, "biogrid.interactions.txt"),
                                 iid = paste0(interactionDataFolder, "iid.interactions.txt"),
                                 mint = paste0(interactionDataFolder, "mint.interactions.txt"),
                                 intact = paste0(interactionDataFolder, "intact.interactions.txt"))

mappingFileLocations <- list(hgnc.map = paste0(interactionDataFolder, "hgnc.mapping.txt"),
                             ensembl.map = paste0(interactionDataFolder, "ensembl.mapping.txt"))

interactionDataList <- list()
for (db in names(interactionFileLocations)) {
  
  if (db == "mint") {
    set.name <- c("uniprot.A", "uniprot.B", "intact.A", "intact.B", "alt.A", "alt.B", "detection.method", "author", "pubmed", "taxid.A", "taxid.B", "interaction.type", "source.database", "interaction.id", "score")
    interactionDataList[[db]] <- read.delim(interactionFileLocations[[db]], stringsAsFactors = FALSE, col.names = set.name)
  }
  interactionDataList[[db]] <- read.delim(interactionFileLocations[[db]], stringsAsFactors = FALSE)
}

mappingFileList <- list()
for (src in names(mappingFileLocations)) {
  mappingFileList[[src]] <- to.char.data.frame(read.delim(mappingFileLocations[[src]], stringsAsFactors = FALSE))
}
colnames(mappingFileList$hgnc.map) <- c("hgnc.symbol", "synonym", "refseq", "entrez.external", "entrez", "uniprotkb.external", "refseq.external", "ensembl.external", "ensembl")
colnames(mappingFileList$ensembl.map) <- c("ensembl", "hgnc.symbol")

mappingFileList$hgnc.map.internal <- mappingFileList$hgnc.map[, c("hgnc.symbol", "synonym", "refseq", "entrez", "uniprotkb.external", "ensembl")]
mappingFileList$hgnc.map.external <- mappingFileList$hgnc.map[, c("hgnc.symbol", "synonym", "refseq.external", "entrez.external", "uniprotkb.external", "ensembl.external")]
mappingFileList$hgnc.map <- NULL


# trim data from downloaded local files

# only use the following id types
used.idtypes <- c("uniprotkb", "ensembl", "refseq", "entrez")

# construct list of id type names
# will not try to find idtype name for these databases 
# as gene symbols were provided
db.with.symbols <- list(biogrid = c(A = 8, B = 9),
                        iid = c(A = 3, B = 4),
                        innatedb = c(A = 5, B = 6))

interactionDataCols <- list(reactome = list(A = c(1,2,3), B = c(4,5,6), interaction.type = 7),
                            innatedb = list(A = c(1,3,5), B = c(2,4,6), symbol.A = 5, symbol.B = 6, detection.method = 7, interaction.type = 12, source.database = 13),
                            intact = list(A = c(1,3), B = c(2,4), detection.method = 7, interaction.type = 12, source.database = 13),
                            biogrid = list(A = c(2,4,6), B = c(3,5,7), symbol.A = 8, symbol.B = 9, source.database = 24),
                            iid = list(A = 1, B = 2, symbol.A = 3, symbol.B = 4, detection.method = 5, source.database = 7, evidence.type = 8, targeting.drugs = 251),
                            mint = list(A = c(1,3), B = c(2,4), detection.method = 7, interaction.type = 12, source.database = 13))

# get processed interaction data frame
idTypeNameList <- list()
newInteractionDataList <- list()
for (db in names(interactionDataList)) {
  
  process.result <- processInteractionData(db)
  idTypeNameList[[db]] <- process.result$id.type
  newInteractionDataList[[db]] <- to.char.data.frame(process.result$new.data.frame)
  rm(process.result)
}

processed.data.folder <- "~/R/"

save.new.int.data(newInteractionDataList, processed.data.folder)
####################################################################################
# need to manually filter out valid data: method of detection, evidence type, etc. #
####################################################################################
# whether database has detection method
# reactome innatedb   intact  biogrid      iid     mint 
# FALSE     TRUE     TRUE    FALSE     TRUE     TRUE 
# man.cur   psi       psi     expr        evid    psi
#           +         +                   *       +
# "+" = psi; "*" = evid

# load methods of detection

prediction.methods <- unlist(strsplit(read_file("prediction.method.name.txt"), "\n"))
prediction.methods <- prediction.methods[prediction.methods != ""]
experimental.methods <- unlist(strsplit(read_file("experimental.method.name.txt"), "\n"))
experimental.methods <- experimental.methods[experimental.methods != ""]

# filter data for detection.method
exprDataList <- list()
for (db in c("innatedb", "intact", "mint")) {
  int.data <- newInteractionDataList[[db]]
  rows.del <- numeric(length = 0)
  rows.keep <- numeric(length = 0)
  for (method in prediction.methods) {
    grep.result <- grep(method, int.data$detection.method)
    rows.del <- c(rows.del, grep.result)
  }
  for (method in experimental.methods) {
    grep.result <- grep(method, int.data$detection.method)
    rows.keep <- c(rows.keep, grep.result)
  }
  cat(db, '---------\n')
  cat("will delete", length(rows.del), "rows\n")
  # print("prediction methods used")
  # print(unique(int.data$detection.method[rows.del]))
  # print("experimental method used")
  # print(unique(int.data$detection.method[rows.keep]))
  print("these extra methods not found")
  print(unique(int.data$detection.method[c(1:length(int.data$detection.method))[-c(rows.del, rows.keep)]]))
  print("assuming the extra methods are experimental")
  exprDataList[[db]] <- to.char.data.frame(int.data[!(c(1:length(int.data[,1])) %in% rows.del), ])
}

# filter detection method for iid
exprDataList$iid <- to.char.data.frame(newInteractionDataList$iid[grep("exp", newInteractionDataList$iid$evidence.type), ])
# reactome is manually curated (used directly)
exprDataList$reactome <- to.char.data.frame(newInteractionDataList$reactome)
# biogrid contains only experimental data (every row has "Experimental.System(.Type)")
exprDataList$biogrid <- to.char.data.frame(newInteractionDataList$biogrid)

for (x in names(exprDataList)) {print(x); print(dim(exprDataList[[x]]))}

# filter for source database
selfSourceName <- c(innatedb = "MI:0974(innatedb)",
                    intact = "psi-mi:MI:0469(IntAct)",
                    mint = "psi-mi:MI:0471(MINT)",
                    biogrid = "BIOGRID",
                    iid = "iid")
# reactome - from self database

selfSourceDataList <- list()
for (db in names(exprDataList)) {
  int.data <- exprDataList[[db]]
  if (db == "iid") {
    use.row <- sapply(strsplit(as.character(int.data$source.database), ";"), FUN = function(x) return("iid" %in% x))
  } else if (db == "reactome") {
    use.row <- c(1:length(int.data[,1]))
  } else if (db %in% names(selfSourceName)) {
    use.row <- which(as.character(int.data$source.database) == selfSourceName[[db]])
  }
  selfSourceDataList[[db]] <- to.char.data.frame(int.data[use.row, ])
}

for (x in names(selfSourceDataList)) {cat(x, "---------\n"); print(dim(selfSourceDataList[[x]]))}

# gene id-symbol conversion
mappingCols <- list(hgnc.map = c(symbol = 1, synonym = 2, refseq = 3, entrez = 4, uniprotkb = 5, ensembl = 6),
                    ensembl.map = c(ensembl = 1, symbol = 2))

# map gene ids
symbolDataList <- list()
for (db in names(selfSourceDataList)) {
  
  int.data <- selfSourceDataList[[db]]
  if (!(db %in% names(db.with.symbols))) {
    for (intr in c("A", "B")) {
      symbol.cols <- symbolForDatabase(int.data[, paste(names(idTypeNameList[[db]]), intr, sep = ".")])
      cn <- colnames(symbol.cols)
      cn[length(cn)] <- paste("symbol", intr, sep = ".")
      colnames(symbol.cols) <- cn
      int.data <- cbind(int.data, symbol.cols)
    }
  }
  
  symbolDataList[[db]] <- to.char.data.frame(int.data)
}

# modify other data bases with symbol info
for (intr in c("A", "B")) {
  
  # extract symbols from innatedb data
  ids <- symbolDataList$innatedb[, paste0("symbol.", intr)]
  
  for (i in 1:length(ids)) {
    hgnc.label <- grep("hgnc", strsplit(ids[i], "\\|")[[1]], value = TRUE)
    
    if (length(hgnc.label) == 0) {
      ids[i] <- NA
    } else {
      ids[i] <- strsplit(hgnc.label, ":|\\(")[[1]][2]
    }
  }
  symbolDataList$innatedb[, paste0("symbol.", intr)] <- ids
  
  # trim biogrid and iid data for symbols existent in hgnc mapping file
  ids <- symbolDataList$biogrid[, paste0("symbol.", intr)]
  ids[!(ids %in% mappingFileList$hgnc.map.internal$hgnc.symbol)] <- NA
  symbolDataList$biogrid[, paste0("symbol.", intr)] <- ids
  
  ids <- symbolDataList$iid[, paste0("symbol.", intr)]
  ids[!(ids %in% mappingFileList$hgnc.map.internal$hgnc.symbol)] <- NA
  symbolDataList$iid[, paste0("symbol.", intr)] <- ids
}

save.new.int.data(symbolDataList, processed.data.folder)

# to remove NA's and empty ""
interactionPairs <- data.frame(symbol.A = character(0), symbol.B = character(0))
for (int.data in symbolDataList) {
  interactionPairs <- distinct(rbind(interactionPairs, int.data[, c("symbol.A", "symbol.B")]))
  interactionPairs <- interactionPairs[(!is.na(interactionPairs[,1])) & (interactionPairs[,1] != ""), ]
  interactionPairs <- interactionPairs[(!is.na(interactionPairs[,2])) & (interactionPairs[,2] != ""), ]
}

# to combine synonyms
syn.count <- 0
multi.syn.list <- list()
pb <- txtProgressBar(min = 0, max = length(interactionPairs[,1]), style = 3)
for (i in 1:length(interactionPairs[,1])) {
  for (j in 1:2) {
    interactionPairs[i, j]
    syn <- unique(mappingFileList$hgnc.map.internal$hgnc.symbol[grep(interactionPairs[i, j], mappingFileList$hgnc.map.internal$synonym)])
    if (length(syn) == 1) {
      interactionPairs[i, j] <- syn
      syn.count <- syn.count + 1
    } else if (length(syn) == 0) {next}
    else {
      # cat("for", interactionPairs[i, j], ":\n")
      # cat("multiple synonyms found:", syn, "will use first\n")
      multi.syn.list[[interactionPairs[i, j]]] <- syn
      interactionPairs[i, j] <- syn[1]
      syn.count <- syn.count + 1
    }
  }
  setTxtProgressBar(pb, value = i)
}
cat(syn.count, "synonyms replaced\n")
interactionPairs <- distinct(interactionPairs)
print(dim(interactionPairs))

# for data saving purpose
int.pairs <- list(int.pair = interactionPairs)
save.new.int.data(int.pairs, processed.data.folder)

# build influnce graph
buildInfluenceGraph <- function(gene.pairs, symmetric = TRUE) {
  if (symmetric) {
    gene.a <- sort(unique(c(gene.pairs[,1], gene.pairs[,2])))
    gene.b <- gene.a
  } else {
    gene.a <- sort(unique(gene.pairs[, 1]))
    gene.b <- sort(unique(gene.pairs[, 2]))
  }
  
  gene.a.name.to.ord <- as.numeric(1:length(gene.a))
  names(gene.a.name.to.ord) <- gene.a
  gene.b.name.to.ord <- as.numeric(1:length(gene.b))
  names(gene.b.name.to.ord) <- gene.b
  cat("building influence graph from gene.a", length(gene.a), "and gene.b", length(gene.b), "\n")
  comm.genes <- intersect(gene.a, gene.b)
  
  #add comm.genes in the use.pairs for a diagonal matrix (A[x,x] = 1)
  if (symmetric) {
    use.pairs <- distinct(data.frame(gene.A = c(gene.pairs[,1], gene.pairs[,2]),
                                     gene.B = c(gene.pairs[,2], gene.pairs[,1]), stringAsFactor = FALSE))
  } else {
    use.pairs <- distinct(data.frame(gene.A = c(gene.pairs[,1]),
                                     gene.B = c(gene.pairs[,2]), stringAsFactor = FALSE))
  }
  
  inf.graph <- Matrix::sparseMatrix(i = gene.a.name.to.ord[use.pairs[,1]], 
                                    j = gene.b.name.to.ord[use.pairs[,2]], 
                                    x = 1, 
                                    dims = c(length(gene.a), length(gene.b)), 
                                    dimnames = list(gene.a, gene.b), 
                                    use.last.ij = TRUE)

  print("preview 1:5, 1:5")
  print(inf.graph[1:5, 1:5])
  # for (i in 1:length(gene.pairs[,1])) inf.graph[gene.pairs[i,1], gene.pairs[i,2]] <- 1
  return(inf.graph)
}

influenceGraph <- as.matrix(buildInfluenceGraph(interactionPairs, symmetric = TRUE))
# for data saving purpose
# influenceGraphBuilt <- list(inf.graph = as.matrix(influenceGraph))
# save.new.int.data(influenceGraphBuilt, processed.data.folder)

# overall data flow: interactionDataList --> newInteractionDataList
# --> exprDataList --> selfSourceDataList

# for (v in ls()) {
#   if (object.size(v) < 100 * 1024 * 1024) save.to.drive <- c(save.to.drive, v)
# }

# DESTIN_FOLDER <- "/Volumes/YUNKAI_HONS/DriverNet/"
# DESTIN_FOLDER <- "~/R/"

# to save
# save.image(paste0(DESTIN_FOLDER, Sys.Date(), "getInput.RData"))

save(list = c("influenceGraph"), file = "out/influence.graph.RData")
ggiOnly <- influenceGraph
save(list = c("ggiOnly"), file = "../GoldDataUseWhenEverPossible/influence.graph.ggi.RData")
########### END OF "generateInfluenceGraph.R" ###################                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            