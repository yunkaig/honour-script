# visualisation of dawnrank results
# plot heatmap with side row annotations on stage/grade of patient
# tile colour = variant/CNV/Mix
library(pheatmap)

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

getEventDF <- function(el) {
  event.df <- lapply(el, function(e) {
    x <- NULL
    for (n in c("cnv", "variant")) {
      if (nrow(e[[n]]) == 0) {next}
      tmp <- e[[n]][, c("GAMuT_ID", "SYMBOL"), drop = F]
      tmp$Event <- n
      rownames(tmp) <- NULL
      x <- rbind(x, tmp)
    }
    return(x)
  })
  event.df <- do.call(rbind, event.df)
  for (i in 1:ncol(event.df)) {event.df[, i] <- as.character(event.df[, i])}
  event.df <- ddply(distinct(event.df), .(GAMuT_ID, SYMBOL), function(d) {
    if (all(c("cnv", "variant") %in% d$Event)) {
      d <- d[1, , drop = F]
      d$Event <- "mixed"
    }
    return(d)
  })
  event.df$Event <- factor(event.df$Event, levels = c("variant", "cnv", "mixed"))
  return(event.df)
}

plotMutation <- function(df, annotation, annot.col, method = "euclidean", 
                       cluster_method = "complete", ...) {
  
  df$event.num <- as.numeric(df$Event)
  for.plot <- dcast(df, GAMuT_ID ~ SYMBOL, value.var = "event.num", fill = NA)
  rownames(for.plot) <- for.plot[, 1]
  for.plot <- as.matrix(for.plot[, -1])
  
  df$code <- ifelse(as.character(df$Event) == 'mixed', 2, 1)
  for.clst <- dcast(df, GAMuT_ID ~ SYMBOL, value.var = "code", fill = -100)
  rownames(for.clst) <- for.clst[, 1]
  for.clst <- as.matrix(for.clst[,-1])
  
  row.clst <- hclust(dist(for.clst, method = method), method = cluster_method)
  col.clst <- hclust(dist(t(for.clst), method = method), method = cluster_method)
  
  # get subset annotations
  annotation <- annotation[row.clst$labels[row.clst$order], , drop = F]
  for.plot <- for.plot[rownames(for.clst), colnames(for.clst)]
  
  p <- pheatmap(for.plot, cluster_rows = row.clst, cluster_cols = col.clst,
                annotation_row = annotation, annotation_colors = annot.col, 
                legend_labels = c("variant", "cnv", "mixed"), 
                legend_breaks = c(1, 2, 3), ...)
  
  
  return(p)
  
}


genes.plot <- genes.plot[1:25]
# el <- getDriverEvent(genes.plot, id.name = "Sample")
el <- getDriverEvent(genes.plot, cnv = valid.cnv, var = valid.var, id.name = "GAMuT_ID")
event.df <- getEventDF(el)

# remove GAMUT prefix 
event.df$GAMuT_ID <- str_replace(event.df$GAMuT_ID, "^GAMuT_", "")
all.results$ID <- str_replace(all.results$Patient, "^GAMuT_", "")

# filter for GAMUT and SYMBOL pairs only
keep <- unlist(apply(event.df, MARGIN = 1, function(e) {
  se1 <- all.results[,"Gene"] == e["SYMBOL"]
  se2 <- all.results[,"ID"] == e["GAMuT_ID"]
  return(all.results$significant[se1 & se2])
}))
event.df <- event.df[as.logical(keep), ]

# get annotation from stage data
stage.data <- read.csv("../All survival_CN_Aug19.csv", as.is = T)
rownames(stage.data) <- stage.data$GAMUT_ID
annot.row <- data.frame(row.names = unique(event.df$GAMuT_ID))
annot.row$Stage <- stage.data[rownames(annot.row), "Stage.group.cat"]
annot.row["844.07G3", "Stage"] <- "I"
annot.row["OV2602.A", "Stage"] <- "III-IV"
annot.row$Stage <- factor(annot.row$Stage, levels = c("NA", "I", "II", "III-IV"))

# remove new EOM case 
eom <- rownames(annot.row)[annot.row$Stage == "EOM"]
annot.row <- annot.row[setdiff(rownames(annot.row), eom), , drop = F]
event.df <- event.df[!(event.df$GAMuT_ID %in% eom), ]

## there is no CA125.high cases
ca125 <- read.csv("../CA125exclude.csv", as.is = T)[,1]
annot.row$CA125.high <- ifelse(rownames(annot.row) %in% ca125, "Yes", "No")
annot.row$CA125.high <- factor(annot.row$CA125.high, levels = c("No", "Yes"))

annot.colour <- list(Stage = c("NA" = "white", "I" = "palegreen", "II" = "yellow", "III-IV" = "tomato"), 
                     CA125.high = c("No" = "grey80", "Yes" = "grey40"))
plotMutation(df = event.df,annotation = annot.row, annot.col = annot.colour, 
             breaks = c(0.5, 1.5, 2.5, 3.5), color = c("steelblue2", "orangered1", "purple"), 
             cellheight = 12, cellwidth = 12, filename = "2019-10-02_hmp_mutplot_stage_allsample.pdf")


