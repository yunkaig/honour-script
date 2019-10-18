library(ggplot2)
library(scales)
library(gridExtra)
library(reshape2)
library(ComplexHeatmap)
library(stringr)
library(ggdendro)
library(pheatmap)

## scripts for analysing GSEA results
# helper functions

getCoreEnrichment <- function(ptw.names, result.dir) {
  return(sapply(ptw.names, FUN = function(e) {
    df <- read.delim(paste0(result.dir, e, ".xls"), sep = "\t", stringsAsFactors = F)
    return(df[df[, "ICHMENT"] == "Yes", "PROBE"])
  }))
}

normalisedGeneCount <- function(genes, all.ptw, subset.ptw, core.genes) {
  # this function gets the gene count in a subset of pathways
  # normalised against k1 and k2 (see below)
  # k1 is the number of pathways that contains a given gene
  k1 <- sapply(searchPathwayList(genes, all.ptw), length)

  # k2 is the number of genes that are in each subset.pathway
  # this assumes each gene in a pathway contributes equally to the function
  # of that pathway
  k2 <- sapply(subset.ptw, length)
  
  gene.search <- searchPathwayList(genes, subset.ptw)
  
  # ngc = normalised.gene.count
  ngc <- sapply(gene.search, FUN = function(e) {
    if (length(e) == 0) return(0)
    return(length(e) * mean(1/k2[e]))
  })
  
  nonzero <- names(k1[k1 != 0])
  ngc[nonzero] <- ngc[nonzero] / k1[nonzero]
  ngc[setdiff(names(ngc), nonzero)] <- 0
  
  if (!is.null(core.genes)) {
    
    # core enriched genes in a large enriched core is penalised
    core.gene.search <- searchPathwayList(genes, core.genes)
    ngc.core <- sapply(core.gene.search, FUN = function(e) {
      if (length(e) == 0) return(0)
      return(length(e) * mean(1/k2[e]))
    })
    
    ngc.core[nonzero] <- ngc.core[nonzero] / k1[nonzero]
    ngc.core[setdiff(names(ngc.core), nonzero)] <- 0
  
  } else {
    ngc.core <- rep(0, length(ngc))
    names(ngc.core) <- names(ngc)
  }
  
  min.val <- min(ngc[ngc != 0])
  ngc.core <- ngc.core / min.val
  ngc <- ngc / min.val - ngc.core[names(ngc)]
  return(list(gene.count = ngc, core.count = ngc.core))
}

getUpsetParts <- function(genes, ptw, report, all.ptw, group.col=NULL, sort.genes=F,
                          low.col = "blue", high.col = "red", export.table = F,
                          table.name = NULL, core.enrich = NULL, dot.size = 2,
                          deg.res = NULL, y.expand = 0.2) {
  # produces parts for combineUpsetParts() for UPset plot
  # export.table = T when want to export conversion table of pathwayNum and pathway 
  if (!is.null(group.col)) {
    ## colour genes by groups
    genes <- unlist(group.col, use.names = F)
    y.col <- rep(names(group.col), sapply(group.col, length))
    names(y.col) <- genes
    
  } else {
    y.col <- rep("black", length(genes))
    names(y.col) <- genes
  }
  
  related.ptws <- searchPathwayList(genes, ptw)
  relate.melt <- melt(related.ptws)
  relate.melt <- relate.melt[, c("L1", "value")]
  colnames(relate.melt) <- c("Gene", "Pathway")
  relate.melt$PathwayNum <- as.numeric(relate.melt$Pathway)
  relate.melt$NES <- report[relate.melt$Pathway, "NES"]
  
  data.df <- dcast(relate.melt, Gene ~ PathwayNum, value.var = "NES", fill = 0, drop = F)
  ## clustering did not give short lines, use default order as before
  # data.df$Gene <- NULL
  # clst <- hclust(dist(data.df))
  # 
  # data.df$Gene <- rownames(data.df)
  plot.df <- melt(data.df, id.var = "Gene", variable.name = "PathwayNum", 
                  value.name = "NES")
  
  # plot.df$Gene <- factor(as.character(plot.df$Gene), levels = clst$labels[clst$order])
  
  plot.df$Gene <- factor(as.character(plot.df$Gene), levels = genes)
  plot.df <- plot.df[order(plot.df$Gene), ]
  plot.df$NES[plot.df$NES == 0] <- NA

  # gene.count for left bar plot
  
  gene.count <- ddply(plot.df[!is.na(plot.df$NES), ], .(Gene), function(e) {
    data.frame(Gene = e[1,"Gene"], meanNES = mean(e$NES))
  })
  ngc.list <- normalisedGeneCount(genes, all.ptw, ptw, core.enrich)
  gene.count[, "count"] <- ngc.list$gene.count[as.character(gene.count[, 'Gene'])]
  gene.count <- ddply(gene.count, .(count), function(e) {
    return(e[order(e$meanNES, decreasing = F), ])
  })
  
  if (sort.genes) {
    # change the order of genes to sort-by-counts
    gene.count$count <- (ngc.list$gene.count + 
                           ngc.list$core.count)[as.character(gene.count[, 'Gene'])]
    gene.count <- gene.count[order(gene.count$count), ]
    genes <- as.character(gene.count$Gene)
    gene.count$Gene <- factor(as.character(gene.count$Gene), levels = genes)
    plot.df$Gene <- factor(as.character(plot.df$Gene), levels = genes)
    gene.count$count <- ngc.list$gene.count[as.character(gene.count[, 'Gene'])]
  }
  
  plot.df <- plot.df[order(plot.df$Gene), ]
  seg.plot <- ddply(plot.df[!is.na(plot.df$NES), ], .(PathwayNum), function(e) {
    return(data.frame(PathwayNum = e[1, "PathwayNum"], 
                      bot = e[1, "Gene"],
                      top = e[nrow(e), "Gene"],
                      NES = e[1, "NES"]))
  })
  
  plot.df$Pathway <- vapply(plot.df$PathwayNum, FUN.VALUE = "1", FUN = function(e) {
    return(as.character(relate.melt[relate.melt$PathwayNum == e, "Pathway"][1]))
  })
  
  deg.df <- gene.count[, c("Gene", "count")]
  deg.df$logFC <- ""
  deg.df$direction <- NA
  if (!is.null(deg.res)) {
    deg.df <- ddply(deg.df, .(Gene), function(e) {
      r <- deg.res[deg.res$Gene == e$Gene, ,drop = F]
      if (nrow(r)  == 0) {
        e$logFC <- "*"
        e$direction <- "exclude"
      } else if (r[, "adj.P.Val"] < 0.05) {
        e$logFC <- round(r[, "logFC"], 2)
        # ignoring the possibility of logFC being 0
        e$direction <- ifelse(e$logFC > 0, "up", "down")
      }
      return(e)
    })
  }
  deg.df$direction <- as.factor(deg.df$direction)
  
  # initialise for core enrichment
  plot.df$CoreEnrich <- F
  gene.count$CoreEnrich <- F
  
  if (!is.null(core.enrich)) {
    for (n in names(core.enrich)) {
      plot.df[(plot.df$Pathway == n) & 
              (plot.df$Gene %in% core.enrich[[n]]), "CoreEnrich"] <- T
    }
  
    core.count <- gene.count
    core.count$count <- ngc.list$core.count[as.character(core.count[, "Gene"])]
    core.count$meanNES <- NA
    core.count$CoreEnrich <- T
    if (any(core.count$Gene != gene.count$Gene)) {
      stop("gene name order inconsistent in core.count and gene.count")
    }

    gene.count <- rbind(core.count, gene.count)
  
  }
  
  
  # ptw.count for top bar plot
  ptw.count <- ddply(plot.df[!is.na(plot.df$NES), ], .(PathwayNum), function(e) {
    p <- relate.melt[relate.melt$PathwayNum == e[1, "PathwayNum"], "Pathway"][1]
    cnt <- nrow(e) / length(ptw[[as.character(p)]])
    ## here the count is normalised by how many genes is present in a pathway
    data.frame(PathwayNum = e[1, "PathwayNum"], count = cnt, meanNES = e[1, "NES"])
  })
  
  ## (deprecated) get pathway from pathway num 
  # seg.plot$Pathway <- vapply(seg.plot$PathwayNum, FUN = function(e) {
  #   return(plot.df[plot.df$PathwayNum == e, "Pathway"][1])
  # }, FUN.VALUE = "A")
  # seg.plot$Pathway <- factor(seg.plot$Pathway, levels = seg.plot$Pathway)
  # 
  # ptw.count$Pathway <- vapply(ptw.count$PathwayNum, FUN = function(e) {
  #   return(plot.df[plot.df$PathwayNum == e, "Pathway"][1])
  # }, FUN.VALUE = "A")
  # ptw.count$Pathway <- factor(ptw.count$Pathway, levels = ptw.count$Pathway)
  # 
  
  ## dp is the main upset plot
  dp <- ggplot() + 
    geom_point(data = plot.df, aes(x = PathwayNum, y = Gene, colour = NES), size = dot.size) +
    geom_point(data = plot.df[plot.df$CoreEnrich, ], aes(x = PathwayNum, y = Gene), 
               shape = 21, colour = "steelblue2", stroke = dot.size/2, size = dot.size) + 
    geom_segment(data = seg.plot, aes(x = PathwayNum, y = top,
                                      xend = PathwayNum, yend = bot, colour = NES)) +
    scale_colour_gradient(na.value = "grey80", low = low.col, high = high.col) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
          panel.grid = element_blank(), 
          axis.text.y = element_text(colour = y.col[as.character(sort(unique(plot.df$Gene)))]))

  ## bp.by.ptw is the bar plot showing how many drivers genes are contained in each pathway
  
  bp.by.ptw <- ggplot() +
    geom_bar(data = ptw.count, aes(x = PathwayNum, y = count, fill = meanNES), 
             stat = "identity", width = 0.8) +
    scale_fill_continuous(na.value = "grey80", low = low.col, high = high.col) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)))
  
  ## bp.by.gene is the bar plot showing how many pathways are associated with each gene
  # print(gene.count)
  
  deg.df$Gene <- factor(deg.df$Gene, levels = levels(gene.count$Gene))
  
  bp.by.gene <- ggplot() + 
    geom_bar(data = gene.count, aes(x = Gene, y = count, fill = meanNES), 
             stat = "identity", position = "stack") + 
    geom_text(data = deg.df, aes(x = Gene, y = count+0.5, label = logFC, 
                                 colour = direction), hjust = 1, size = 3) +
    scale_fill_gradient(na.value = "grey60", low = low.col, high = high.col) + 
    scale_colour_manual(values = c(up = "red", down = "blue", exclude = "green")) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    coord_flip() + scale_y_reverse(expand = expand_scale(mult = c(y.expand, 0)))
  
  if (export.table) {
    if (is.null(table.name)) {
      table.name <- paste0(Sys.Date(), "_PathwayNum.txt")
    }
    write.table(distinct(relate.melt[order(relate.melt$PathwayNum), c("PathwayNum", "Pathway")]), file = table.name, 
                col.names = T, row.names = F, sep = "\t", quote = F)
  }
  
  return(list(upset = dp, 
              bp.ptw = bp.by.ptw,
              bp.gene = bp.by.gene,
              gene.count = gene.count, 
              deg.df = deg.df))
  
}

combineUpsetParts <- function(parts, main="upset", left="bp.gene", top="bp.ptw", 
                              pad = 2, top.adj = pad + 16, 
                              title = "UpSet plot") {
  # combine the parts created by getUpsetParts() function and produce the upset plot
  main <- parts[[main]]
  left <- parts[[left]]
  top <- parts[[top]]
  
  # get legend from main
  tmp <- ggplot_gtable(ggplot_build(main))
  n <- sapply(tmp$grobs, function(e) e$name)
  main.lgd <- tmp$grobs[[which(n == "guide-box")]]
  
  ## cleaning 
  main.clean <- main +
    theme(legend.position = "none", axis.title.y = element_blank(), 
          axis.text.y = element_text(hjust = 0.5), axis.ticks.y = element_blank(),
          plot.margin = unit(c(pad, 4, 4, pad), "pt"))
  
  top.clean <- top + 
    theme(legend.position = "none", axis.title.x = element_blank(), 
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          plot.margin = unit(c(4, 4, pad, top.adj), "pt"), 
          panel.background = element_blank())
  
  left.clean <- left + 
    theme(legend.position = 'none', axis.title.y = element_blank(), 
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
          plot.margin = unit(c(pad, pad, 4, 4), "pt"), panel.background = element_blank())
  
  return(grid.arrange(main.lgd, top.clean, left.clean, main.clean, ncol = 2, nrow = 2,
               widths = unit(c(100, 650), "pt"), heights = unit(c(150, 300), "pt"),
               top = title))
}

plotUpset <- function(genes, ptw, report, all.ptw, group.col = NULL, sort.genes = F, 
                      core.enrich = NULL, dot.size = 2, pad = 2, top.adj = NULL, 
                      title = paste(deparse(substitute(genes)), "UpSet plot"), 
                      filename = NULL, deg.res = NULL, y.expand = 0.2) {
  if (is.null(filename)) {
    filename <- paste(Sys.Date(), deparse(substitute(genes)), "Upset.pdf", sep = "_")
  }
  table.name <- paste0(Sys.Date(), "_PathwayNum_", deparse(substitute(genes)), ".txt")
  parts <- getUpsetParts(genes = genes, ptw = ptw, report = report, all.ptw = all.ptw, 
                         group.col = group.col, sort.genes = sort.genes, 
                         export.table = T, table.name = table.name, 
                         core.enrich = core.enrich, dot.size = dot.size, 
                         deg.res = deg.res, y.expand = y.expand)
  
  # top.adj + 2 is to compensate for PDF plotting 
  # (somehow the margin changes a bit when using pdf)
  if (is.null(top.adj)) {
    top.adj <- (max(str_length(genes)) - 1) * 4
  }
  p <- combineUpsetParts(parts, pad = pad, top.adj = top.adj + 2, title = title)
  pdf(file = filename, width = 12, height = 9, onefile = T)
  grid.arrange(p)
  dev.off(dev.list()["pdf"])
  cat("file saved to", filename, "\n")
}

countFreq <- function(core, all.ptw) {
  count.df <- count(melt(core), "value")
  count.df$percent <- count.df$freq / length(core)
  count.df$value <- as.character(count.df$value)
  count.df <- count.df[order(count.df$freq, decreasing = T), ]
  rownames(count.df) <- NULL
  return(count.df)
}

pathwayToMatrix <- function(ptw) {
  all.g <- unique(unlist(ptw))
  emp.df <- rep(0, length(all.g))
  names(emp.df) <- all.g
  ptw.mat <- sapply(ptw, function(p) {
    d <- emp.df
    d[p] <- 1
    return(d)
  })
  return(t(ptw.mat))
}

clusterPathway <- function(ptw, method = "binary", min_freq = 0, y.expand= 0.1, 
                           cluster.genes = T) {
  
  ptw.mat <- pathwayToMatrix(ptw)
  clst <- hclust(dist(ptw.mat, method = method))
  dendr <- dendro_data(clst)
  # filter for genes occurred in more than one pathways
  ptw.mat <- ptw.mat[, colSums(ptw.mat) >= min_freq]
  
  # dg <- ggplot() + 
  #   geom_segment(data = dendr$segments, aes(x=x,y=y,xend=xend,yend=yend,
  #                                           colour=colour)) + 
  #   scale_y_reverse(expand = expand_scale(add = c(0.01, y.expand))) + 
  #   geom_text(data = dendr$labels, aes(x=x,y=y-0.1, colour=colour, label = label), 
  #             hjust=0, size = 2.5) +
  #   scale_colour_discrete(na.value = "black") +
  #   scale_x_continuous(expand = expand_scale(add = c(0.1, 0.1))) + coord_flip() +
  #   theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
  #         axis.ticks = element_blank(), axis.text.x = element_blank(), 
  #         axis.text.y = element_blank(), legend.position = "none",
  #         panel.background = element_rect(fill = NA), panel.grid = element_blank())
  
  ptw.melt <- melt(ptw.mat, varnames = c("pathway", "gene"))
  ptw.melt$pathway <- factor(ptw.melt$pathway, levels = clst$labels[clst$order])
  if (cluster.genes) {
    clst.gene <- hclust(dist(t(ptw.mat), method = method))
    ptw.melt$gene <- factor(ptw.melt$gene, levels = clst.gene$labels[clst.gene$order])
  }
  
  rast <- ggplot() + geom_raster(data = ptw.melt, aes(x = gene, y = pathway, 
                                                      fill = as.factor(value))) + 
    scale_fill_manual(values = c("0" = "grey40", "1" = "white")) + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position = "none")
  # plot <- grid.arrange(dg, rast, ncol = 2, nrow = 1, widths = unit(c(300, 800), "pt"), 
  #                      heights = unit(800, "pt"))
    
  return(rast)
}


#####
## !! deprecated because takes too much resource to run
#####
dendr.as.list <- function(dendr, top.y = NULL) {
  # this function converts dendr.object to a nested list recursively
  # e.g. list(1, list(2, 3, list(4, 5)))
  # each time, the functions returns list(A,B, dendr.as.list(C.top.y))
  
  # note all(dendr$segments$y > dendr$segments$yend) == TRUE

  if (is.null(top.y)) {
    vert.seg <- dendr$segments[dendr$segments$x == dendr$segments$xend, ]
    hori.seg <- dendr$segments[dendr$segments$y == dendr$segments$yend, ]
    top.y <- max(hori.seg$y)
  } else {
    chs <- dendr$segments$y <= top.y
    vert.seg <- dendr$segments[dendr$segments$x == dendr$segments$xend & chs, ]
    hori.seg <- dendr$segments[dendr$segments$y == dendr$segments$yend & chs, ]
  }

  stem <- lapply(top.y, function(e) {
    tmp <- hori.seg[hori.seg$y == e, ]
    # assumming xend contains the extreme values of a bar
    return(c(y = tmp[, "y"][1], 
             xmax = max(tmp[, "xend"]),
             xmin = min(tmp[, "xend"])))
  })
  
  all.res <- list()
  for (st in stem) {
    res <- list()
    # for each stem, find the vertical bars "ends" that is 
    # directly under that stem
    chs1 <- vert.seg$y == st["y"]
    chs2 <- vert.seg$x <= st["xmax"]
    chs3 <- vert.seg$x >= st["xmin"]
    
    # new vertical bars
    tmp.2 <- vert.seg[chs1 & chs2 & chs3, ]
    
    chs4 <- tmp.2$yend == 0
    # get ends, yend == 0
    res <- c(res, tmp.2[chs4 , "x"])
    
    # get others to make the recursive call
    for (y in tmp.2[!chs4, "yend"]) {
      res <- c(res, list(dendr.as.list(dendr, top.y = y)))
    }
    
    all.res <- c(all.res, res)
  }
  
  return(all.res)
  
}

groupByLevel <- function(seg, labels, lv = 1) {

  group <- c()
  
  hori.seg <- seg[seg$y == seg$yend & seg$x != seg$xend, ]
  vert.seg <- seg[seg$x == seg$xend & seg$y != seg$yend, ]
  
  # get immediate vertical bars
  top.y <- max(hori.seg$y)
  tmp <- vert.seg[vert.seg$y == top.y, ]
  if (lv == 0) {
    # find minimum xend
    vb.min <- tmp[tmp$x == min(tmp$x), ]
    while (vb.min$yend != 0) {
      next.hb <- hori.seg[hori.seg$x == vb.min$xend & hori.seg$y == vb.min$yend, ]
      vb.min <- vert.seg[vert.seg$x == min(next.hb$xend) & 
                           vert.seg$y == next.hb$yend[1], ]
    }
    xmin <- vb.min$xend
    vb.max <- tmp[tmp$x == max(tmp$x), ]
    while (vb.max$yend != 0) {
      next.hb <- hori.seg[hori.seg$x == vb.max$xend & hori.seg$y == vb.max$yend, ]
      vb.max <- vert.seg[vert.seg$x == max(next.hb$xend) & 
                           vert.seg$y == next.hb$yend[1], ]
    }
    xmax <- vb.max$xend
    
    lb <- as.character(labels[labels$x <= xmax & labels$x >= xmin, "label"])
    # lb <- as.character(xmin:xmax)
    group[lb] <- 1
    return(group)
  }
  
  x <- tmp[tmp$yend == 0, "x"]
  lb <- as.character(labels[labels$x %in% x, "label"])
  # lb <- as.character(x)

  ## adding "leftouts" to the group
  group[lb] <- 1
  
  # get the possible next y values
  next.y <- tmp$yend
  next.y <- next.y[next.y != 0]
  
  for (y in next.y) {
    chs1 <- seg$y < y
    chs2 <- !(seg$x %in% x & seg$yend == 0)
    new.seg <- seg[chs1 & chs2, ]
    group <- c(group, groupByLevel(new.seg, labels, lv - 1) + max(group))
  }
  n <- as.character(sort(factor(names(group), levels = labels$label)))
  return(group[n])
  
}

# msigdb report ------
# folder where index.html sits
result.folder <- "./run.separate/all_c2_functional/out/Jul24/my_analysis.Gsea.1563942369049/"
pheno <- "MOC"
report.file <- dir(result.folder, pattern = paste0("gsea_report_for_", pheno, "_[[:digit:]].*xls"))
report <- read.delim(paste0(result.folder, report.file), stringsAsFactors = F, sep = '\t')
ptw.file <- paste0(result.folder, "edb/gene_sets.gmt")
ptw <- parsePathway(ptw.file)
names.clean <- sapply(strsplit(names(ptw), "\\|"), FUN = function(e) e[[1]])
names(ptw) <- names.clean       # the number of significant genes

# filter with FDR < 0.25, nom-p < 0.05
ori.report <- report
report <- report[(report$NOM.p.val <= 0.05) & (report$FDR.q.val <= 0.25), ]
rownames(report) <- as.character(report$NAME)
nrow(report)
sig.ptw <- ptw[report$NAME]
# flat.sig.ptw <- unlist(sig.ptw, use.names = F)
# inc.genes <- unique(flat.sig.ptw)
# freq.count <- sapply(inc.genes, FUN = function(e) length(which(flat.sig.ptw == e)))

# load("../GitHubDriverNet/predicted.drivers_JULY11.RData")
# com.genes <- intersect(genes.plot, dawn.genes)
genes.plot <- read.table(text="1    TP53
2    KRAS
                         3    CFTR
                         4   APLP2
                         5     MYC
                         6   LRRK2
                         7  CDKN2A
                         8   ERBB2
                         9    FLNA
                         10  CDK12
                         11   MDM2
                         12 RECQL4
                         13  ARIH1
                         14 TGFBR1
                         15  SMAD3
                         16    FGB
                         17 TRIM41
                         18     AR
                         19  IKZF3
                         20 NEDD4L
                         21   NPM1
                         ", stringsAsFactors = F)[,2]
dawn.genes <- read.table(text="1    TP53
2    KRAS
                         3  CDKN2A
                         4   ERBB2
                         5   CDK12
                         6     JUP
                         7   PSMD3
                         8  PIK3CA
                         9    ACLY
                         10   RARA
                         11  RPL19
                         12    TTN
                         13  MED24
                         14  LRRK2
                         15  DDX3X
                         16  MED12
                         17  RPL23
                         18    MYC
                         19  ERBB3
                         20  SUZ12", stringsAsFactors = F)[,2]

uni.col <- list(tomato = intersect(dawn.genes, genes.plot), 
                green = setdiff(genes.plot, dawn.genes),
                blue = setdiff(dawn.genes, genes.plot))

plotUpset(genes = uni.genes, ptw = sig.ptw, report = report, group.col = uni.col, 
          sort.genes = T)

uni.top.col <- list(tomato = top.genes, palegreen2 = uni.genes)
plotUpset(genes = c(uni.genes, top.genes), ptw = sig.ptw, all.ptw = ptw, report = report, 
          group.col = uni.top.col,sort.genes = T, core.enrich = core.enrich)

# plot UpSet plot for common results between singscore and GSEA
load("../singscore/2019-10-13_comm.ptw.RData")
plotUpset(genes = c(uni.genes, top.genes), ptw = sig.ptw[unlist(comm.ptw, use.names = F)], 
          all.ptw = ptw, report = report, group.col = uni.top.col,sort.genes = T, 
          core.enrich = core.enrich[unlist(comm.ptw, use.names = F)], top.adj = 18, 
          filename = "2019-10-13_uni.genes_top.genes_comm.ptw.pdf")

load("../singscore/2019-10-13_sig.phe.RData")
# test if GSEA and singscore results independent and intersection due to chance
N <- length(ptw) # 3520
n.singscore <- length(unique(unlist(sig.phe[2:4]))) # 778
n.gsea <- length(sig.ptw)

n.gp.sp <- length(unlist(comm.ptw)) # gsea positive, singscore positive 33
n.gp.sn <- n.gsea - n.gp.sp # 123
n.gn.sp <- n.singscore - n.gp.sp # 745
n.gn.sn <- N - n.gp.sp - n.gp.sn - n.gn.sp # 2619

# use chi-sq test
test.mat <- matrix(c(n.gp.sp, n.gn.sp, n.gp.sn, n.gn.sn), 2, 2,
                   dimnames = list(c("GSEA+", "GSEA-"), c("singscore+", "singscore-")))
test.res <- chisq.test(as.table(test.mat))
save(list = c("test.mat", "test.res"), file = "2019-10-13_chisq-test-gsea-singscore.RData")
############# main ends --------

# set report to reactome.unique ------


# folder where index.html sits
result.folder <- "./outJun25.uni.react/my_analysis.Gsea.1561434648118/"
pheno <- "MOC"
report.file <- dir(result.folder, pattern = paste0("gsea_report_for_", pheno, "_[[:digit:]].*xls"))
report <- read.delim(paste0(result.folder, report.file), stringsAsFactors = F, sep = '\t')
ptw.file <- paste0(result.folder, "edb/gene_sets.gmt")
ptw <- parsePathway(ptw.file)
names.clean <- sapply(strsplit(names(ptw), "\\|"), FUN = function(e) e[[1]])
names(ptw) <- names.clean

# filter with FDR < 0.25, nom-p < 0.05
report <- report[(report$NOM.p.val <= 0.05) & (report$FDR.q.val <= 0.25), ]
rownames(report) <- as.character(report$NAME)
nrow(report)
sig.ptw <- ptw[report$NAME]


load("../GitHubDriverNet/predicted.drivers_JULY11.RData")
com.genes <- intersect(genes.plot, dawn.genes)
uni.genes <- union(genes.plot, dawn.genes)
related.ptws <- searchPathwayList(com.genes, sig.ptw)
relate.melt <- melt.list(related.ptws)
relate.melt <- relate.melt[, c("L1", "value")]
colnames(relate.melt) <- c("Gene", "Pathway")
relate.melt$PathwayNum <- as.numeric(relate.melt$Pathway)
relate.melt$NES <- report[relate.melt$Pathway, "NES"]
top.plot <- levels(relate.melt$Pathway)
plot.df <- relate.melt[relate.melt$Pathway %in% top.plot, ]


dp <- ggplot() + 
  geom_point(data = plot.df, aes(x = PathwayNum, y = Gene)) +
  # geom_point(data = relate.melt, aes(x = PathwayNum, y = NES), colour = "tomato") + 
  # scale_y_continuous(limits = c(min(relate.melt$NES), NA), oob = rescale_none) + 
  theme(axis.text.x = element_text(angle = 0))
# dp

bp <- ggplot() + 
  geom_point(data = plot.df, aes(x = PathwayNum, y = NES), colour = "tomato")
# bp
grid.arrange(dp, bp, nrow = 2, ncol = 1)

nes.vs.ptw <- ggplot(data = cbind(ori.report[, c("NES", "NOM.p.val", "FDR.q.val")], ptwn = 1:nrow(ori.report))) + 
  geom_point(aes(x = ptwn, y = NES), 
             colour = ifelse(ori.report$NOM.p.val <= 0.05 & ori.report$FDR.q.val <= 0.25, "red", "grey40"), 
             position = "identity", stat = "identity") + theme(legend.position = "none")
nes.vs.ptw

comp.lst <- list()
for (p in list(react, msig)) {
  temp <- c(length = length(p))
  temp <- c(temp, vapply(com.genes, FUN.VALUE = 0L, FUN = function(e) {
    length(searchPathwayList(e, p)[[1]])
  }))
  comp.lst <- c(comp.lst, list(temp))
}
names(comp.lst) <- c('Reactome', "MSigDB")
comp.df <- as.data.frame(comp.lst)
comp.df

comp.df$stat <- rownames(comp.df)
comp.plot <- data.frame(stat = rep(comp.df$stat, 2), count = c(comp.df$Reactome, comp.df$MSigDB), 
                        db = c(rep("Reactome", nrow(comp.df)), rep("MSigDB", nrow(comp.df))), 
                        length = c(rep(comp.df["length", "Reactome"], nrow(comp.df)), rep(comp.df["length", "MSigDB"], nrow(comp.df))),
                        stringsAsFactors = F)
comp.plot$stat <- factor(comp.plot$stat, levels = comp.df$stat)
comp.plot <- comp.plot[order(comp.plot$stat), ]

comp.bar <- ggplot(data = comp.plot[-c(1,2), ]) + 
  geom_bar(aes(x = stat, y = count/length * 100 , fill = factor(db)), 
           stat = "identity", position = "dodge") + 
  scale_fill_discrete(name = "Database") + ylab("normalised count") +
  theme(axis.text.x = element_text(angle = 90))
comp.bar


# gene count for highly enriched pathways ------


gene.count <- numeric()
for (g in unlist(sig.ptw, use.names = F)) {
  if (g %in% names(gene.count)) {
    gene.count[[g]] <- gene.count[[g]] + 1
  }
  else {
    gene.count[[g]] <- 1
  }
}
# this line to account for when the gene to plot are not present in the significant pathways
gene.count[uni.genes[!(uni.genes %in% names(gene.count))]] <- 0

count.df <- data.frame(Gene = names(gene.count), Count = gene.count)
count.df <- count.df[order(count.df$Count, decreasing = T), ]
count.df$Gene <- factor(count.df$Gene, levels = count.df$Gene)
count.df$PR <- percent_rank(nrow(count.df):1)
count.df$tool <- "None"
count.df[count.df$Gene %in% genes.plot, 'tool'] <- "DriverNet"
count.df[count.df$Gene %in% dawn.genes, 'tool'] <- "DawnRank"
count.df[count.df$Gene %in% com.genes, 'tool'] <- "Both"
count.df$tool <- factor(count.df$tool, levels = c("Both", "DriverNet", "DawnRank", "None"))
count.df$seq <- 1:nrow(count.df)


plot.data <- count.df[uni.genes, ]
plot.data <- plot.data[order(plot.data[, "Count"]), ]
plot.data$Gene <- factor(as.character(plot.data$Gene), levels = as.character(plot.data$Gene))

# plot the numeber of pathways containing each driver gene
bp <- ggplot(data = plot.data) + 
  geom_point(aes(x = Count, y = Gene, colour = tool)) + 
  geom_text(aes(x = Count, y = Gene, label = Count), hjust = -0.5, angle = 0) + 
  theme(axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5,
                                   hjust = 1))
bp

# plot the common pathways between driver genes
ptw.list <- searchPathwayList(uni.genes, sig.ptw)
overlap.list <- lapply(ptw.list, FUN = function(p1) {
  sapply(ptw.list, FUN = function(p2) {
    if (length(p1) == 0 | length(p2) == 0) return(0)
    length(which(p1 %in% p2))
  })
})
overlap.df <- as.data.frame(overlap.list)
for(x in 1:nrow(overlap.df)) overlap.df[x, x] <- 0
overlap.df[1:4,1:4]

# linkage between driver genes and jess genes
gene.comb <- lapply(uni.genes, function(e) return(data.frame(uni = e, top = top.genes)))
gene.comb <- do.call("rbind",gene.comb)
gene.comb[,1] <- as.character(gene.comb[,1])
gene.comb[,2] <- as.character(gene.comb[,2])
n <- apply(gene.comb, MARGIN = 1, FUN = function(e) {
  s <- sapply(sig.ptw, FUN = function(x) {
    return((e[1] %in% x) & (e[2] %in% x))
  })
  nms <- sort(names(sig.ptw[s]))
  return(data.frame(count = length(nms),
           names = paste(nms, collapse = ";")))
})

gene.comb <- cbind(gene.comb, do.call("rbind", n))
for (i in 1:ncol(gene.comb)) {
  gene.comb[, i] <- as.character(gene.comb[, i])
}
gene.comb <- gene.comb[order(gene.comb$count, decreasing = T), ]
write.table(gene.comb, "2019-09-02_gene_overlap.txt", sep = "\t", col.names = T, 
            row.names = F, quote = F)

x <- unique(unlist(strsplit(gene.comb$names, split = ";")))
len <- sapply(x, function(e) length(ptw[[e]]))
len <- sort(len)
write.table(data.frame(length=len, name = names(len)), 
            "2019-09-02_sort_overlap_pathway.txt", 
            col.names = T, row.names = F, quote = F, sep = "\t")

ptw.num.df <- data.frame(pathway = sort(names(sig.ptw)), 
                         pathway.num = 1:length(sig.ptw), 
                         stringsAsFactors = F)

ptw.cluster <- parsePathway("ptw.cluster.txt", F, " ")
folder <- "ptw.cluster.plots/"
for (n in names(ptw.cluster)) {
  file1 = paste0(folder, n, '_cluster_genes.pdf')
  file2 = paste0(folder, n, "_increasing_freq.pdf")
  mat <- pathwayToMatrix(sig.ptw[ptw.cluster[[n]]])
  if (ncol(mat) > 200) {
    s <- colSums(mat) > 1
    mat <- mat[,s]
  }
  
  nr <- nrow(mat)
  nc <- ncol(mat)
  w <- nc/6+7
  h <- nr/5+5
  pheatmap(mat, clustering_distance_cols = "binary", clustering_distance_rows = "binary",
           color = c("palegreen2", "tomato"), cellwidth = 8, cellheight = 8, fontsize = 8,
           height = h, width = w, filename = file1, silent = T)
  ord <- order(colSums(mat))
  pheatmap(mat[,ord], clustering_distance_cols = "binary", 
           cluster_cols = F, color = c("palegreen2", "tomato"), 
           cellwidth = 8, cellheight = 8, fontsize = 8,
           height = h, width = w, filename = file2, silent = T)
  
}

get.closest.ptw <- function(name, dist.mat = dm, n = 5) {
  return(sort(dist.mat[name, ])[2:(n+1)])
}

get.closest.ptw("SMID_BREAST_CANCER_LUMINAL_A_DN")
get.closest.ptw("OUYANG_PROSTATE_CANCER_MARKERS")

# the threshold distance

thres <- 0.7
# determine the closest pathways to each pathway

close.ptws <- apply(dm, MARGIN = 1, function(e) {
  return(names(e)[e > 0 & e < thres])
})
close.ptws <- close.ptws[sapply(close.ptws, length) > 0]






