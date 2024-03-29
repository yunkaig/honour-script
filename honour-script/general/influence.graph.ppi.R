for (pkg in c("igraph",
              "tidyverse", 
              "stringr",
              "Matrix")) {
  if (!requireNamespace(pkg)) {
    install.packages(pkg, repos = "https://cran.ms.unimelb.edu.au/", quiet = TRUE)
  }
  library(pkg, character.only = TRUE)
}


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

buildInfluenceGraph <- function(gene.pairs, symmetric) {
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



## read in mapping files
interactionDataFolder <- "interactionData/"
mappingFileLocations <- list(hgnc.map = paste0(interactionDataFolder, "hgnc.mapping.txt"),
                             ensembl.map = paste0(interactionDataFolder, "ensembl.mapping.txt"))
mappingFileList <- list()
for (src in names(mappingFileLocations)) {
  mappingFileList[[src]] <- to.char.data.frame(read.delim(mappingFileLocations[[src]], stringsAsFactors = FALSE))
}
colnames(mappingFileList$hgnc.map) <- c("hgnc.symbol", "synonym", "refseq", "entrez.external", "entrez", "uniprotkb.external", "refseq.external", "ensembl.external", "ensembl")
colnames(mappingFileList$ensembl.map) <- c("ensembl", "hgnc.symbol")

## read in ppi
print("Reading ppi from file ------")
ppi.folder <- "mocNe-yg-sl/"
ppi.files <- list(iid = paste0(ppi.folder, "graph_IID_lcs.graphml"),
                  huri = paste0(ppi.folder, "graph_HuRI_III_lcs.graphml"))
                  # string = paste0(ppi.folder, "graph_STRING_lcs.graphml"))

ppi.list <- list()
for (db in names(ppi.files)) {
  ppi.list[[db]] <- igraph::read.graph(ppi.files[[db]], format = "graphml")
}

# get pairwise interactions (long_data_frame)
print("Running paired.long.df -----")
paired.long.df <- list()
for (db in names(ppi.list)) {
  cat("working on", db, "\n")
  paired.long.df[[db]] <- igraph::as_long_data_frame(ppi.list[[db]])[c("from_name", "to_name")]
  cat(db, "completed\n")
}

# convert to symbols
print("Running symbols.long.df -----")
symbols.long.df <- list()
for (db in names(paired.long.df)) {
  
  # Ignoring any number after '-', e.g. Q9H2S6-1 will be converted to Q9H2S6
  # I guess that these dashes might mean protein transcripts, 
  # but still have to confirm with Sam
  ldf <- paired.long.df[[db]]
  ldf[,1] <- str_remove(ldf[,1], "-[:digit:]$")
  ldf[,2] <- str_remove(ldf[,2], "-[:digit:]$")
  
  cat("working on", db, "\n")
  conv.df <- mappingFileList$hgnc.map[, c("uniprotkb.external", "hgnc.symbol")]
  conv.df <- conv.df[which(conv.df[,1] != ""), ]
  symb.ldf <- ldf
  for (i in 1:2) {
    symb.ldf[,i] <- geneSymbolFromID(convDF = conv.df, ids = symb.ldf[,i])
  }
  symb.ldf <- na.omit(symb.ldf)
  
  # to combine synonyms
  syn.count <- 0
  multi.syn.list <- list()
  
  for (i in 1:nrow(symb.ldf)) {
    for (j in 1:2) {
      symb.ldf[i, j]
      syn <- unique(mappingFileList$hgnc.map.internal$hgnc.symbol[grep(symb.ldf[i, j], mappingFileList$hgnc.map.internal$synonym)])
      if (length(syn) == 1) {
        symb.ldf[i, j] <- syn
        syn.count <- syn.count + 1
      } else if (length(syn) == 0) {next}
      else {
        # cat("for", symb.ldf[i, j], ":\n")
        # cat("multiple synonyms found:", syn, "will use first\n")
        multi.syn.list[[symb.ldf[i, j]]] <- syn
        symb.ldf[i, j] <- syn[1]
        syn.count <- syn.count + 1
      }
    }
    
  }
  cat("\n", syn.count, "synonyms replaced\n")
  symb.ldf <- distinct(na.omit(symb.ldf))
  
  symbols.long.df[[db]] <- symb.ldf
  cat(db, "completed\n")
}

# combine ppi's

combined.ppi <- data.frame(from_name = character(0), to_name = character(0))
for (db in names(symbols.long.df)) {
  combined.ppi <- rbind(combined.ppi, symbols.long.df[[db]])
}
combined.ppi <- distinct(combined.ppi)

influenceGraph <- buildInfluenceGraph(combined.ppi, symmetric = TRUE)
ppiOnly <- influenceGraph
save.image(file = paste0("out/", Sys.Date(), "_symbol_long_data_frame.RData"))
save(list = "ppiOnly", file = "../GoldDataUseWhenEverPossible/influence.graph.ppi.RData")
                                                                                                                                                                                                                                                                                                                                                                                                                                        