# functions to parse the reactome pathway mapping file
# generate a large list
library(stringr)
library(readr)

parsePathway <- function(file, info.col = T, sep = "\t") {
  ptw.text <- read_file(file)
  ptw.text <- strsplit(ptw.text, "\n")[[1]]
  ptw.text <- strsplit(ptw.text, sep)
  
  ptw.list <- lapply(ptw.text, FUN = function(e) {
    return(e[3: length(e)])
  })
  if (info.col) {
    ptw.names <- lapply(ptw.text, FUN = function(e) {
      return(paste(e[1:2], collapse = "|"))
    })
  } else {
    ptw.names <- lapply(ptw.text, FUN = function(e) {
      return(paste(e[1]))
    })
  }
  names(ptw.list) <- ptw.names
  return(ptw.list)
}

writePathway <- function(ptw, file, info.col = TRUE) {
  text <- sapply(ptw, FUN = function(p) {
    return(paste0(p, collapse = "\t"))
  })
  if (info.col) {
    name.cols <- sub("\\|", "\t", names(ptw))
  } else {
    name.cols <- names(ptw)
  }
  
  text <- paste(name.cols, text, sep = "\t")
  write_lines(text, file)
}

to.char.df <- function(df) {
  for (i in 1:ncol(df)) {
    df[, i] <- as.character(df[, i])
    
  }
  return(df)
}

searchPathwayList <- function(genes, ptw) {
  # returns a list(gene = ptw.name)
  res.list <- lapply(genes, FUN = function(g) {
    res.v <- sapply(ptw, FUN = function(e) {
      g %in% e
    })
    names(res.v)[res.v]
  })
  names(res.list) <- genes
  return(res.list)
}

is.all.sorted <- function(ptw) {
  !all(sapply(ptw, FUN = function(e) {is.unsorted(e)}))
}

findCommonRoots <- function(symbols, single.root = FALSE, use.names = FALSE) {
  # single.root set to TRUE when force to find ONE common root
  # use.names set to TRUE to get vector with symbols as names and 
  # the corresponding root as values
  # this function is used to predict common roots of a vector of strings (gene symbols) 
  
  if (single.root) {
    if (use.names) {
      warning("use.names argument will be set to FALSE as single.root is set to TRUE")
    }
    # use brute force to find a single root
    root = character(0)
    for (i in 1:min(str_length(symbols))) {
      i.char <- str_sub(symbols, i, i)
      if (all(i.char == i.char[1])) {
        root <- c(root, i.char[1])
      } else {
        return(paste0(root, collapse = ""))
      }
    }
  } else {
    # remove digits at the end
    no.end.digits <- str_remove(symbols, "\\d{1,}$")
    if (!use.names) return(unique(no.end.digits))
    names(no.end.digits) <- symbols
    return(no.end.digits)
  }
}

replaceWithRoot <- function(symbols, root.list = NULL, show.names = FALSE) {
  # show.names is TRUE only for debug uses
  if (length(root.list) == 0) {
    root.list <- findCommonRoots(symbols, use.names = TRUE)
  }
  root.list <- c(root.list, "")


  # search in the symbols from left to right
  # when encounters a symbols with a different root
  # replace the previous symbols with their common root
  # will not replace when only one subtype occurs
  # e.g. ERBB2 stays ERBB2 because no other ERBB's before or after
  root.count <- 0     # counts how many symbols with the same root occurs consecutively
  root.len <- length(root.list) - 1
  replaced.symbols <- character(0)
  next.root <- root.list[1]
  
  for (i in 1:root.len) {
    root <- next.root
    next.root <- root.list[i + 1]
    root.count <- root.count + 1
    
    # if the next root differs (and counter > 1), then append replaced.symbols and reset counter
    if (root != next.root) {
      if (root.count > 1) {
        replaced.symbols <- c(replaced.symbols, root)
        
      } else {
        elem <- names(root)
        names(elem) <- elem
        replaced.symbols <- c(replaced.symbols, elem)
        
      }
      root.count <- 0
    }
  }
  
  if (!show.names) names(replaced.symbols) <- NULL
  return(replaced.symbols)
}


getSingleIndependentPathway <- function(ptw) {
  # expand independent set from the first pathway
  hub.name <- names(ptw)[1]
  hub <- ptw[[1]]
  taken <- hub.name
  p.names <- names(ptw)
  last.hub <- hub
  new.genes <- hub
  while (length(new.genes) > 0) {
    p.names <- setdiff(names(ptw), taken)
    for (n in p.names) {
      p <- ptw[[n]]
      if (any(p %in% new.genes)) {
        taken <- c(taken, n)
        hub.name <- c(hub.name, n)
        hub <- union(hub, p)
        diff <- TRUE
      }
    }
    new.genes <- setdiff(hub, last.hub)
    last.hub <- hub

  }
  return(list(names = hub.name, genes = hub))
}

getIndependentPathways <- function(ptw, check.valid = FALSE) {
  # returns a list of independent gene sets 
  # each element in ind.ptw --> two subelements: $pathway.names & $genes
  # self.check for debugging use: to check if overlapping exists between gene sets
  # ind.ptw used to store pathway names only
  ind.ptw <- list()
  while (length(ptw) > 0) {
    ind.ptw[[length(ind.ptw) + 1]] <- getSingleIndependentPathway(ptw)
    ptw <- ptw[setdiff(names(ptw), ind.ptw[[length(ind.ptw)]]$names)]
  }
  
  if (check.valid && length(ind.ptw) > 1) {
    # check for overlap between indepdent gene sets

    for (i in 1:(length(ind.ptw) - 1)) {
      for (j in (i + 1):length(ind.ptw)) {
        if (any(ind.ptw[[i]]$genes %in% ind.ptw[[j]]$genes)) {
          warning(paste("found overlapping genes between independent sets", i, j))
        }
        if (any(ind.ptw[[i]]$names %in% ind.ptw[[j]]$names)) {
          warning(paste("found overlapping names between independent sets", i, j))
        }
      }
    }
  }
  
  return(ind.ptw)
}

duplicatedPathways <- function(ptw) {
  dup.ptw <- duplicated(ptw)
  dup.names <- names(ptw)[dup.ptw]
  cat("Found", length(dup.names), "duplicated pathways in", deparse(substitute(ptw)), "\n")
  return(dup.names)
}

rmIncomplete <- function(ptw, panel.genes) {
  i <- sapply(ptw, function(e) {
    all(e %in% panel.genes)
  })
  cat("removed:", length(which(!i)), "pathways\n")
  cat(length(which(i)), "remains\n")
  return(ptw[i])
}

# gate the trimmed gene sets at minimum length of 15 (number from GSEA.out file)
trimIncomplete <- function(ptw, panel.genes, min.len = 15) {
  ptw <- lapply(ptw, "intersect", panel.genes)
  return(ptw[sapply(ptw, "length") >= min.len])
}



