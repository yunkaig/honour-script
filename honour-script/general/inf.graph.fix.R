library(Matrix.utils)
library(dplyr)

## fix the influence graph matrix for the "|" (pipe) in the gene names issue

cast.sparse.mat <- function(df, row.var, col.var, value.var) {
  if (mode(df[, row.var]) == "character") {
    df[, row.var] <- as.factor(df[, row.var])
    row.lvs <- levels(df[, row.var])
  } else {
    row.lvs <- 1:max(df[, row.var])
  }
  if (mode(df[, col.var]) == "character") {
    df[, col.var] <- as.factor(df[, col.var])
    col.lvs <- levels(df[, col.var])
  } else {
    col.lvs <- 1:max(df[, col.var])
  }
  
  i <- as.numeric(df[, row.var])
  j <- as.numeric(df[, col.var])
  mat <- sparseMatrix(i = i, j = j, x = df[, value.var],
                      dimnames = list(row.lvs, col.lvs))
  return(mat)
}

fix.pipe <- function(x) {
  summ <- summary(x)
  ln.df <- data.frame(A = rownames(x)[summ$i], 
                      B = colnames(x)[summ$j], 
                      value = 1, stringsAsFactors = F)
  rownames(ln.df) <- NULL
  for (g1 in c("A", "B")) {
    for (i in c("A", "B")) {
      ln.df[, i] <- as.character(ln.df[, i])
    }
    g2 <- ifelse(g1 == "A", "B", "A")
    se <- grep("\\|", ln.df[, g1])
    # extract rows with pipe in g1 column
    tmp <- ln.df[se, ]
    ln.df <- ln.df[-se, ]
    
    tmp.lst <- apply(tmp, MARGIN = 1, function(e) {
      e <- unlist(e)
      gsep <- strsplit(e[g1], "\\|")[[1]]
      b <- e[g2]
      names(b) <- NULL
      res <- data.frame(g1 = gsep, g2 = b, value = 1)
      colnames(res) <- c(g1, g2, "value")
      return(res)
      
    })
    tmp.df <- do.call(rbind, tmp.lst)
    rownames(tmp.df) <- NULL
    ln.df <- distinct(rbind(ln.df, tmp.df))
  }
  
  return(cast.sparse.mat(ln.df, "A", "B", "value"))
  
}


load("../GoldDataUseWhenEverPossible/influence.graph.all.RData")

ggiOnly <- fix.pipe(ggiOnly)
ppiOnly <- fix.pipe(ppiOnly)
ggiANDppi <- fix.pipe(ggiANDppi)

save(list = c("ggiOnly", "ppiOnly", "ggiANDppi"), 
     file = "../GoldDataUseWhenEverPossible/influence.graph.all.RData")

