library(reshape2)
library(ggplot2)
library(plyr)
library(stringr)
library(gridExtra)
library(ComplexHeatmap)
library(ggdendro)
library(tidyverse)


# helper functions
getDriverEvent <- function(genes, cnv = mocCNVData, var = mocVariantData, 
                           comb.cnv = T, comb.var = F, id.name = "GAMuT_ID") {
  
  sym.split <- strsplit(as.character(cnv$Gene.Symbols), split = ", ")
  # made a mistake here to use grep(g, mocCNVData$Gene.Symbols) to account for the 
  # comman separation, should split and use exact equal "==" 
  # This BUG resulted in the high "AR" CNV count when there is actually only 2 CNVs
  if (id.name != "GAMuT_ID") {
    cnv.cn <- colnames(cnv)
    var.cn <- colnames(var)
    cnv <- cnv[, setdiff(cnv.cn, "GAMuT_ID")]
    var <- var[, setdiff(var.cn, "GAMuT_ID")]
    cnv.cn <- colnames(cnv)
    var.cn <- colnames(var)
    cnv.cn[cnv.cn == id.name] <- "GAMuT_ID"
    var.cn[var.cn == id.name] <- "GAMuT_ID"
    colnames(cnv) <- cnv.cn
    colnames(var) <- var.cn
  }
  el <- list()
  for (g in genes) {
    i <- sapply(sym.split, FUN = function(e) g %in% e)
    cnv.info <- cnv[i, c("GAMuT_ID", "Chromosome", "Event")]
    if (nrow(cnv.info) > 0) {
      cnv.info$SYMBOL <- g
      cnv.info$Event <- factor(cnv.info$Event, 
                               levels = c("High Copy Gain", "Homozygous Copy Loss"))
      levels(cnv.info$Event) <- c("HA", "HD")
      if (comb.cnv) {
        # combine HA and HD mix in CNV
        cnv.info <- ddply(distinct(cnv.info), .(GAMuT_ID), function(e) {
          if (nrow(e) == 1) return(e)
          e <- e[1, ]
          e$Event <- "Mix"
          return(e)
        })
      }
    }
    el[[g]] <- list(cnv = cnv.info,
                    variant = var[var$SYMBOL == g, 
                                  c("GAMuT_ID", "CHROM", "Consequence_Rank", "SYMBOL")])
    if (comb.var) {
      # combine consequence rank mix in variant
      el[[g]]$variant <- ddply(distinct(el[[g]]$variant), .(GAMuT_ID), function(e) {
        if (nrow(e) == 1) return(e)
        r <- paste0(sort(e[, "Consequence_Rank"]), collapse = ",")
        e <- e[1, ]
        e$Consequence_Rank <- r
        return(e)
      })
    }
  }
  el
}

countEvents <- function(g, de = eventList[[g]]) {
  # de = driver.event = eventList[[g]] (The event for THIS gene)
  cnv.event <- integer()
  for (e in de$cnv$Event) {
    if (!(e %in% names(cnv.event))) cnv.event[e] <- 1
    else cnv.event[e] <- cnv.event[e] + 1
  }
  var.cons <- integer()
  for (r in de$variant$Consequence_Rank) {
    r <- as.character(r)
    if (!(r %in% names(var.cons))) var.cons[r] <- 1
    else var.cons[r] <- var.cons[r] + 1
  }
  list(cnv = cnv.event, variant = var.cons)
}

event.string <- function(event.count) {
  cnv.event <- event.count$cnv
  var.cons <- event.count$variant
  paste0("[", paste(names(cnv.event), cnv.event, sep = ":", collapse = ","), "][",
         paste(paste0("variant", names(var.cons)), var.cons, sep = ":", collapse = ","), "]")
}

getAdjGene <- function(g, inf.graph = influenceGraph) {
  if (!isSymmetric.matrix(inf.graph)) {
    stop("inf.graph argument must be symmetrical")
  }
  if (!(g %in% rownames(inf.graph))) {
    stop(paste("genes not found in inf.graph:", g[!(g %in% rownames(inf.graph))]))
  }
  return(colnames(inf.graph)[inf.graph[g,] == 1])
}

plotEventCountBar <- function(genes = names(eventList), el = eventList, lab.size = 3, 
                              gene.lab = F, horizontal = F) {
  # el = eventList -- originally used eventCount, but one patient can have multiple events in one gene
  # remove genes with no info in el
  if (any(!(genes %in% names(el)))) {
    warning(paste(paste(genes[!(genes %in% names(el))], collapse = ", "), 
                  "are not found in eventList, will be ignored"))
    genes <- genes[genes %in% names(el)]
  }
  # plot bar graph regardless of subtypes of variant or CNV
  # patients with both varian and CNV will be represented by overlapping of stacked columns
  # plot.mat e.g TP53 Variant = 16 means "there are 16 samples with variant ONLY, excluding overlaps with CNV"
  plot.mat <- vapply(genes, FUN.VALUE = c(Variant = 1, CNV = 1, Overlap = 1),
                     FUN = function(e) {
                       var <- unique(el[[e]]$variant$GAMuT_ID)
                       cnv <- unique(el[[e]]$cnv$GAMuT_ID)
                       overlap <- intersect(var, cnv)
                       var <- setdiff(var, overlap)
                       cnv <- setdiff(cnv, overlap)
                       return(c(Variant = length(var),
                                CNV = length(cnv), Overlap = length(overlap)))
                     })
  plot.mat <- as.data.frame(t(plot.mat), stringsAsFactors = FALSE)
  plot.mat$Gene <- rownames(plot.mat)
  melt.mat <- melt(plot.mat, id.vars = "Gene")
  melt.mat <- melt.mat[melt.mat[, "value"] > 0, ]
  
  melt.mat$variable <- factor(melt.mat$variable, levels = c("CNV", "Overlap", "Variant"))
  melt.mat$Gene <- factor(melt.mat$Gene, levels = genes)
  # add genes with zero counts
  if (any(!(genes %in% melt.mat$Gene))) {
    melt.mat <- rbind(melt.mat, data.frame(Gene = genes[!(genes %in% melt.mat$Gene)],
                                           value = 0, variable = "CNV"))
  }
  
  melt.mat <- ddply(melt.mat, .(Gene), function(e) {
    e <- e[order(e$variable), ]
    e$lab.y = c(0, cumsum(e$value[nrow(e):1]))[nrow(e):1] + e$value/2
    return(e)
  })
  
  col.plot <- c("tomato", rgb(255, 200, 129, maxColorValue = 255), "palegreen2")
  txt.ang <- ifelse(horizontal, 0, 90)
  bp <- ggplot() + 
    geom_bar(data = melt.mat, aes(x = Gene, y = value, fill = variable), stat = "identity") + 
    geom_text(data = melt.mat, aes(x = Gene, y = lab.y, label = value), 
              size = lab.size, angle = txt.ang) + 
    scale_fill_manual(name = "Event", values = col.plot) + 
    ylab("Number of Samples") + xlab("Predicted Driver Genes") + 
    theme(axis.text.x = element_text(angle = txt.ang, vjust = 0.5, hjust = 1), 
          axis.ticks.x = element_blank())
  
  
  if (gene.lab){
    gene.lab.mat <- ddply(melt.mat, .(Gene), function(e) {
      return(data.frame(Gene = e[1, "Gene"], gene.y = sum(e$value), stringsAsFactors = F))
    })
    bp <- bp + 
      geom_text(data = gene.lab.mat, aes(x = Gene, y = gene.y, label = Gene), 
                hjust = -.1, size = lab.size, angle = txt.ang) +
      scale_y_continuous(expand = expand_scale(mult = c(0, 0.15)))
  }
  
  if (horizontal) {return(bp + coord_flip())}
  return(bp)
  
}


# Deprecated, see differential expression R file with limma-voom method
# plotDifferentialExpression <- function(genes, 
#                                        tumour = t(mocTumourExpression), 
#                                        normal = t(benExpression)) {
#   # produce a violin plot of the differential expression of genes
#   # for driverNet, need to t() tumour and normal as input
#   # because need patients to be colnames
#   
#   if (any(!(genes %in% rownames(tumour)))) {
#     warning(paste("Some genes not found in expression matrices, ignoring:",
#             paste(genes[!(genes %in% rownames(tumour))], collapse = ", ")))
#     genes <- genes[genes %in% rownames(tumour)]
#   }
#   tumour <- tumour[genes, ]
#   normal <- normal[genes, ]
#   plot.mat <- tumour - rowMeans(normal)
#   t.ncol <- ncol(tumour)
#   c.ncol <- ncol(normal) + t.ncol
#   signif <- apply(cbind(tumour, normal), MARGIN = 1, FUN = function(e) {
#     t.test(e[1:t.ncol], e[(t.ncol+1):c.ncol])$p.value
#   })
#   signif[signif < 0.001] <- "***"
#   signif[0.001 <= signif & signif < 0.01] <- "**"
#   signif[0.01 <= signif & signif < 0.05] <- "*"
#   signif[signif >= 0.05] <- ""
#   
#   signif <- data.frame(Gene = names(signif), 
#                        value = signif, 
#                        height = apply(plot.mat, 1, max))
#   
#   melt.mat <- melt(plot.mat, varnames = c("Gene", "Patient"))
#   vp <- ggplot() + 
#     geom_hline(yintercept = 0, colour = "grey50") + 
#     geom_violin(data = melt.mat, scale = "width", trim = TRUE, adjust = 1.5, 
#                 draw_quantiles = c(0.25, 0.5, 0.75),
#                 aes(x = Gene, y = value, fill = factor(Gene)), width = 0.6) + 
#     geom_jitter(height = 0, width = 0.1) + 
#     geom_text(data = signif, aes(x = Gene, y = height + 1.5, 
#                                  label = value), 
#               size = 6, angle = 90, vjust = 0.75) + 
#     theme(legend.position = "none", 
#           axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 1),
#           axis.ticks.x = element_blank())
#   vp
# }


plotCNVBar <- function(genes, el = getDriverEvent(genes, cnv = valid.cnv, var = NULL, comb.cnv = T), 
                       lab.size = 3) {
  mat.plot <- do.call("rbind", lapply(el[genes], function(e) e$cnv))
  mat.plot <- ddply(mat.plot, .(SYMBOL, Event), function(e) {
    e <- cbind(e[1, c("SYMBOL", "Event")], count = nrow(e))
    return(e)
  })
  mat.plot <- ddply(mat.plot, .(SYMBOL), function(e) {
    tmp <- e$count
    tmp <- tmp[length(tmp):1]
    tmp <- tmp / 2 + c(0, cumsum(tmp)[-length(tmp)])
    e$lab.y <- tmp[length(tmp):1]
    return(e)
  })
  
  bp <- ggplot(data = mat.plot, aes(x = SYMBOL)) + 
    geom_bar(aes(y = count, fill = Event), stat = "identity",
             position = "stack") + 
    geom_text(aes(y = lab.y, label = count), size = lab.size) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  bp
}

clusterData <- function(df, by = "mutation") {
  
  for.clst <- distinct(df[, c("GAMuT_ID", "SYMBOL", by)])
  dcast.mat <- as.matrix(dcast(for.clst, GAMuT_ID ~ SYMBOL, value.var = by))
  dcast.mat[is.na(dcast.mat)] <- ifelse(by == "mutation", 0, "NONE")
  rownames(dcast.mat) <- as.character(dcast.mat[, "GAMuT_ID"])
  dcast.mat <- dcast.mat[, -1]
  
  if (by == "variant") {
    # for variant clustering, distance calculated by same/not the same
    # assumes the distances (difference) is the same between any two consequences and 
    # wild types
    clst <- hclust(dist2(dcast.mat, pairwise_fun = function(x, y) {
      sum(as.numeric(x != y) ^ 2)
    }))
  } else if (by == "cnv") {
    clst <- hclust(dist2(dcast.mat, pairwise_fun = function(x, y) {
      d <- 0
      for (i in 1:length(x)) {
        a <- x[i]
        b <- y[i]
        if (a == b) next
        if (all(c("NONE", "Mix") %in% c(a, b))) {
          d <- d + 2^2
          next
        }
        if (all(c("HD", "HA") %in% c(a, b))) {
          d <- d + 2^2
          next
        }
        else d <- d + 1
      }
      return(sqrt(d))
    }))
  } else {
    # in cnv and mutation clustering, distance calculated by vector (euclidean)
    clst <- hclust(dist(dcast.mat, method = "euclidean"))
  }
  return(clst)
}

getPlotMatrix <- function(genes, el) {
  el <- el[genes]
  
  varcount <- do.call("rbind", lapply(el, FUN = function(e) e$variant))[, c("GAMuT_ID", "CHROM", "Consequence_Rank", "SYMBOL")]
  # need to change colname CHROM -> Chromosome
  colnames(varcount) <- c("GAMuT_ID", "Chromosome", "Consequence_Rank", "SYMBOL")
  cnvcount <- do.call("rbind", lapply(el, FUN = function(e) e$cnv))
  
  # combine var and cnv count data.frames into mat.plot
  
  mat.plot <- rbind(cbind(varcount, Event = NA), cbind(cnvcount, Consequence_Rank = NA))
  mat.plot <- ddply(distinct(mat.plot), .(GAMuT_ID, SYMBOL), function(e) {
    if (length(which(is.na(e$Consequence_Rank))) == 0) {
      return(e)
    }
    if (length(which(is.na(e$Event))) == 0) {
      return(e)
    }
    tmp <- e[is.na(e$Event), ]
    tmp$Event <- e[!is.na(e$Event), "Event"]
    
    return(tmp)
  })
  
  mat.plot[is.na(mat.plot$Consequence_Rank), "Consequence_Rank"] <- "NONE"
  mat.plot[is.na(mat.plot$Event), "Event"] <- "NONE"
  mat.plot <- distinct(mat.plot)
  # determine mutation status in mat.plot
  cn <- colnames(mat.plot)
  cn[cn == "Event"] <- "cnv"
  cn[cn == "Consequence_Rank"] <- "variant"
  colnames(mat.plot) <- cn
  mat.plot$mutation <- as.numeric(!(mat.plot$variant == "NONE" & mat.plot$cnv == "NONE"))
  mat.plot <- ddply(mat.plot, .(GAMuT_ID), function(e) {
    fill.mat <- data.frame(GAMuT_ID = e[1, "GAMuT_ID"], SYMBOL = genes[!(genes %in% e$SYMBOL)], 
                           Chromosome = "NONE", cnv = "NONE", variant = "NONE", mutation = 0, 
                           stringsAsFactors = F)
    return(rbind(e, fill.mat))
  })
  mat.plot$cnv <- factor(mat.plot$cnv, levels = c("NONE", "HA", "HD", "Mix"))
  var.levels <- unique(mat.plot$variant)
  var.levels <- var.levels[order(nchar(var.levels), var.levels)]
  var.levels <- c("NONE", var.levels[var.levels != "NONE"])
  mat.plot$variant <- factor(mat.plot$variant, levels = var.levels)
  return(mat.plot)
}

plotGenomeAbberationHeatmap <- function(genes, el, cluster.by = "mutation", 
                                        show.cnv = T, dendrogram = F, cut = 3, cols = NULL,
                                        one.piece = F, dendr.top.margin = 6,
                                        dendr.bot.margin = 70, y.expand = 3,
                                        x.expand =  0.01) {
  # plot heatmap for variant and cnv data
  # "split-cells" for multiple types in single individual
  if (any(!(genes %in% names(el)))) {
    warning(paste(c("Some genes not found in expression matrices, ignoring:",
                    genes[!(genes %in% names(el))]), collapse = " "))
    genes <- genes[genes %in% rownames(tumour)]
  }
  if (!(cluster.by %in% c('mutation', "cnv", "variant"))) {
    stop("cluster.by should be one of: 'mutation', 'cnv', 'variant'")
  }
  mat.plot <- getPlotMatrix(genes, el)
  
  clst <- clusterData(df = mat.plot, by = cluster.by)
  cat("clustered sample ids:\n")
  cat('"', paste0(clst$labels[clst$order], collapse = '","'), '" \n', sep = "")
  group <- cutree(clst, k = cut)
  
  h.fac <- 0.7     # this is the height argument in aes for non-split cells
  w.fac <- 0.7
  
  mat.plot$GAMuT_ID <- factor(mat.plot$GAMuT_ID, levels = clst$labels[clst$order])
  mat.plot$SYMBOL <- factor(mat.plot$SYMBOL, levels = genes)
  
  mat.plot <- ddply(mat.plot, .(GAMuT_ID, SYMBOL), function(e) {
    cbind(e, y = as.numeric(e$GAMuT_ID) - 0.5 * h.fac + 
            (0.5 + (nrow(e):1 - 1)) * h.fac/nrow(e),
          h = h.fac/nrow(e))
  })
  
  hmp <- ggplot(data = mat.plot) + 
    geom_tile(aes(x = SYMBOL, y = GAMuT_ID), fill = NA) +    # this is to show the sample IDs on the y-axis
    geom_tile(aes(x = SYMBOL, y = y, 
                  fill = variant, 
                  colour = if (show.cnv) {
                    factor(cnv, levels = c("HA", "HD", "Mix", "NONE"))
                  } else NA,
                  width = w.fac, height = ifelse(h == 1, h.fac * h, h)),
              size = 0.7) + 
    scale_fill_manual(name = "Variant\nConsequence\nRank", 
                      values = c("1" = "red",
                                 "2" = "pink",
                                 "3" = "purple", 
                                 "4" = "orange", 
                                 "NONE" = "grey80"), 
                      labels = c("1" = "truncating frame-shift", 
                                 "2" = "essential splice site", 
                                 "3" = "in-frame indel", 
                                 "4" = "stop codon loss", 
                                 "NONE" = "None"),
                      guide = guide_legend(order = 1)) + 
    scale_colour_manual(name = "CNV Events",
                        values = c("HA" = "green3", 
                                   "HD" = "steelblue1", 
                                   "Mix" = "yellow",
                                   "NONE" = NA), 
                        guide = guide_legend(order = 2, 
                                             override.aes = list(fill = "grey80"))) + 
    ylab("Samples") + xlab("Predicted driver genes") +
    theme(axis.text = element_text(size = 8), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.ticks = element_blank(), panel.grid = element_blank(), 
          axis.text.y = element_text(colour = group[levels(mat.plot$GAMuT_ID)] + 1))
  if (!dendrogram) {
    return(hmp)
  }
  
  dendr <- dendro_data(clst, type = "rectangle")
  # dendr plots horizontally initially
  dendr <- colourDendrogram(dendr, group)
  dg <- ggplot() + 
    geom_segment(data = dendr$segments, aes(x=x,y=y,xend=xend,yend=yend,
                                            colour=colour)) + 
    scale_y_reverse(expand = expand_scale(add = c(0.01, y.expand))) + 
    geom_text(data = dendr$labels, aes(x=x,y=y-0.1,colour=colour, label = label), 
              hjust=0, size = 2.5) +
    scale_colour_discrete(na.value = "black") +
    scale_x_continuous(expand = expand_scale(add = c(x.expand,x.expand))) + coord_flip() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.ticks = element_blank(), axis.text.x = element_blank(), 
          axis.text.y = element_blank(), legend.position = "none",
          panel.background = element_rect(fill = NA), panel.grid = element_blank())
  if (!one.piece) {
    return(list(heatmap = hmp, dendrogram = dg))
  }
  # arrange the plots in one piece
  hmp <- hmp + 
    theme(plot.margin = unit(c(0,4,0,2), "pt") ,
          axis.title.y = element_blank(),
          axis.text.y = element_blank())
  dg <- dg + 
    theme(plot.margin = unit(c(dendr.top.margin, 2, dendr.bot.margin,4), "pt"))
  return(grid.arrange(dg, hmp, ncol = 2, nrow = 1, widths = unit(c(100, 300), "pt")))
  
}

colourDendrogram <- function(dendr, group) {
  # this function returns colourised dendr by group
  cut <- length(unique(group))
  heights <- unique(dendr$segments$y)
  heights <- sort(heights, decreasing = T)
  
  maxh <- mean(heights[cut], heights[cut-1])
  # > all(dendr$segments$yend <= dendr$segments$y)
  # [1] TRUE
  coloured <- dendr$segments$yend <= maxh
  dendr$segments$colour <- NA
  dendr$labels$colour <- group[as.character(dendr$labels$label)]
  
  for (n in unique(group)) {
    g <- group[group == n]
    cat("===> group", n,"<=== \n")
    cat(names(g), sep = '","')
    cat("\n")
    x.val <- dendr$labels$x[dendr$labels$label %in% names(g)]
    min.x <- min(x.val)
    max.x <- max(x.val)
    keep1 <- min.x <= dendr$segments$x & dendr$segments$x <= max.x
    keep2 <- min.x <= dendr$segments$xend & dendr$segments$xend <= max.x
    dendr$segments$colour[keep1 & keep2 & coloured] <- n
  }
  dendr$segments$colour <- as.factor(dendr$segments$colour)
  dendr$labels$colour <- as.factor(dendr$labels$colour)
  return(dendr)
}

plotHeatmapAndBar <- function(genes, el = eventList, cluster.by = "mutation", 
                              tile.size = 4, bp = NULL, gene.lab = F,
                              group.col = NULL, cols = NULL, dendr.top.margin = 2, 
                              dendr.bot.margin = NULL, cut = 3) {
  # group.col in the form list(colour1 = c("gene1", "gene2"...) ...)
  # get sample size
  n.sample <- length(unique(unlist(lapply(el, function(e) {
    return(unique(c(e$cnv$GAMuT_ID, e$variant$GAMuT_ID)))
  }))))
  if (!is.null(group.col)) {
    x.col <- c()
    for (col in names(group.col)) {
      col.tmp <- rep(col, length(group.col[[col]]))
      names(col.tmp) <- group.col[[col]]
      x.col <- c(x.col, col.tmp)
    }
    genes <- names(x.col)
  }
  if (is.null(dendr.bot.margin)) {
    dendr.bot.margin <- 36 + max(str_length(genes)) * 4
  }
  # tile.size is the size of a single tile in the heatmap, in the unit of mm
  if (is.null(bp)) {
    bp <- plotEventCountBar(genes = genes, el = el, lab.size = tile.size*0.5, gene.lab = gene.lab) + 
      theme(legend.key.size = unit(tile.size, "mm"))
  }
  
  hmp.res <- plotGenomeAbberationHeatmap(genes = genes, el = el, cluster.by = cluster.by,
                                         cols = cols, dendr.top.margin = 0, 
                                         dendr.bot.margin = 0, 
                                         one.piece = F, dendrogram = T, cut = cut,
                                         y.expand = 2.5, x.expand = 0.5)
  
  hmp <- hmp.res$heatmap + 
    theme(legend.key.size = unit(tile.size, "mm"), axis.text.y = element_blank())
  if (!is.null(group.col)) {
    hmp <- hmp + theme(axis.text.x = element_text(colour = x.col[genes]))
  }
  
  dg <- hmp.res$dendrogram + 
    theme(plot.margin = unit(c(dendr.top.margin, 0, dendr.bot.margin, 0), "pt"))
  
  hmp.clean <- hmp + 
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
          axis.title.y = element_blank(), legend.position = "none", 
          plot.margin = unit(c(0,0,12,0), "pt"))    # 12 pt bottom for axis.text.x
  tmp <- ggplot_gtable(ggplot_build(hmp))
  leg <- which(sapply(tmp$grobs, function(e) e$name) == "guide-box")
  hmp.lgd <- tmp$grobs[[leg]]
  
  bp.clean <- bp +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
          axis.title.y = element_blank(), axis.text.x = element_blank(), 
          axis.title.x = element_blank(), legend.position = "none", 
          plot.margin = unit(c(0, 0, 0, 0), "pt"), panel.background = element_blank(), 
          panel.grid = element_blank())
  tmp <- ggplot_gtable(ggplot_build(bp))
  leg <- which(sapply(tmp$grobs, function(e) e$name) == "guide-box")
  bp.lgd <- tmp$grobs[[leg]]
  hmp.w <- tile.size * length(genes)
  hmp.h <- tile.size * n.sample
  
  grid.arrange(grob(), bp.clean, bp.lgd, dg, hmp.clean, hmp.lgd, ncol = 3, nrow = 2, 
               heights = unit(c(tile.size*15, hmp.h), "mm"), 
               widths = unit(c(50, hmp.w, tile.size*10), "mm"),
               top = textGrob("Genomic abberations"), padding = unit(5, "mm"))
}

getGrobSize <- function(g) {
  return(c(sum(convertWidth(g$widths, "in", valueOnly = T)),
           sum(convertHeight(g$heights, "in", valueOnly = T))))
}

plotPDF <- function(plot, file = NULL, size = NULL, portrait = TRUE) {
  if (is.null(file)) {
    file <- paste0(Sys.Date(), "_", deparse(substitute(plot)), ".pdf")
  }
  file <- str_replace_all(file, "\'|\"", "")
  # create portrait or landscape
  if (is.null(size)) {
    if ("grob" %in% class(plot)) {
      size <- round(getGrobSize(plot)) + 1
    } else {
      if (portrait) size <- c(7, 9)
      else size <- c(9, 7)
    }
  }
  cat("set PDF size to", size, "\n")
  pdf(file = file, onefile = FALSE, height = size[2], width = size[1])
  
  # use grid.arrange, otherise PDF will be empty
  grid.arrange(plot)
  
  dev.off()
  cat("Plot saved to", file, "\n")
}

check.grobs <- function(plot) {
  # this function iterates through all GRaphic OBjects (grob) of the input plot
  # the plot argument must be a ggplot object
  # all non-empty grobs are displayed separately
  grobs <- ggplot_gtable(ggplot_build(plot))$grobs
  for (i in 1:length(grobs)) {
    x <- grobs[[i]]
    if (x$name == "NULL") next
    grid.arrange(x)
    readline(prompt = paste(i, x$name, "Press enter to continue ..."))
  }
}