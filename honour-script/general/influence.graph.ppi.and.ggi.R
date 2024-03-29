library(Matrix)
library(tidyverse)

adj_matrix_to_long_data_frame <- function(m) {
  ldf <- data.frame(row = character(0), col = character(0))
  
  # determine if to use dimnames 
  if ((NA %in% colnames(m)) || (NA %in% rownames(m))) {
    warning("NA found in colnames/rownames of matrix, may affects result")
  }
  rs <- if (length(rownames(m)) == 0) c(1:nrow(m)) else rownames(m)
  cs <- if (length(colnames(m)) == 0) c(1:ncol(m)) else colnames(m)
  for (e in rs) {
    ones <- cs[m[e,] == 1]
    if (length(ones > 0)) {
      ldf <- rbind(ldf, data.frame(row = rep(e, length(ones)), col = ones))
      
    }
  }
  ldf
}

adj_matrix_to_adj_list <- function(m, margin = 1) {
  # margin = 1 : rownames as keys and colnames as elements
  # margin = 2 : colnames as keys and rownames as elements
  if (margin == 2) {
    m <- t(m)
  } else if (margin != 1) {
    stop("margin argument must be 1 (rows as key) or 2 (cols as key)")
  }
  # determine if to use dimnames 
  if ((NA %in% colnames(m)) || (NA %in% rownames(m))) {
    warning("NA found in colnames/rownames of matrix, may affects result")
  }
  rs <- if (length(rownames(m)) == 0) c(1:nrow(m)) else rownames(m)
  cs <- if (length(colnames(m)) == 0) c(1:ncol(m)) else colnames(m)
  al <- list()
  for (e in rs) {
    ones <- cs[m[e,] == 1]
    if (length(ones > 0)) al[[e]] <- ones
  }
  al
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

to.char.data.frame <- function(df) {
  for (i in 1:length(df[1, ])) df[, i] <- as.character(df[, i])
  return(df)
}

combineInfluenceGraph <- function(inf.1, inf.2) {
  
  inf.1.gene <- colnames(inf.1)
  inf.2.gene <- colnames(inf.2)
  
  comm.genes <- intersect(inf.1.gene, inf.2.gene)
  inf.1.only <- setdiff(inf.1.gene, comm.genes)
  inf.2.only <- setdiff(inf.2.gene, comm.genes)
  all.genes <- sort(union(inf.1.gene, inf.2.gene))
  
  print(data.frame(set = c(paste0(deparse(substitute(inf.1)), "inf.1"), paste0(deparse(substitute(inf.2)), "inf.2"), "comm.genes", "inf.1.only", "inf.2.only", "all.genes"),
                   Ngenes = c(length(inf.1.gene), length(inf.2.gene), length(comm.genes), length(inf.1.only), length(inf.2.only), length(all.genes))))

  
  is.equal <- (inf.1[comm.genes, comm.genes] == inf.2[comm.genes, comm.genes])
  
  # build combined influence graph - inf.comb
  inf.comb <- sparseMatrix(i = c(), j = c(), x = 1, dims = rep(length(all.genes), 2))
  dimnames(inf.comb) <- list(all.genes, all.genes)
  inf.comb[inf.1.gene, inf.1.gene] <- inf.1
  inf.comb[inf.2.gene, inf.2.gene] <- inf.2
  inf.comb
}



# load ggi.inf.graph
load("out/Run on 2-may/influence.graph.RData")
ggi.inf.graph <- influenceGraph.symmetric
rm(list = c("influenceGraph", "influenceGraph.symmetric"))

# load ppi.inf.graph
load("out/13 May 3019 PPI-run-IID-HuRI/influence.graph.RData")
ppi.inf.graph <- influenceGraph
rm(influenceGraph)


# calculate the difference between ggi and ppi
# gene-wise:
commonGenes <- intersect(colnames(ppi.inf.graph), colnames(ggi.inf.graph))
ggiOnly <- setdiff(colnames(ggi.inf.graph), commonGenes)
ppiOnly <- setdiff(colnames(ppi.inf.graph), commonGenes)
cat("\t/PPI:", ncol(ppi.inf.graph), "\\\t\t/GGI:", ncol(ggi.inf.graph), "\\\nPPI only:", length(ppiOnly), "\tCommon:", length(commonGenes), "\tGGI only:", length(ggiOnly), "\n", sep = "")

# combined the two influencegraphs
influenceGraph <- combineInfluenceGraph(ggi.inf.graph, ppi.inf.graph)
save(list = c("influenceGraph", "commonGenes", "ggiOnly", "ppiOnly"), file = "out/influence.graph.RData")

                                                                                                                                                                                                                                                                                                                                                                                                                             