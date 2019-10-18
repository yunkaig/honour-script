# post singscore analysis
# visualisation 
library(gplots)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(plyr)
library(stringr)
library(pheatmap)
library(ggdendro)
library(ggpubr)

# # unfiltered
# dir.file <- "out/2019-07-31_unfiltered_directional_phenotype_scores_all.txt"
# undir.file <- "out/2019-07-31_unfiltered_undir_scores_all.txt"
# dir.scoredf <- read.delim(dir.file, sep = "\t", stringsAsFactors = F)
# undir.scoredf <- read.delim(undir.file, sep = "\t", stringsAsFactors = F)

file.list <- paste0("out/", dir('out/', "^2019-08-05.*_results\\.txt"))
df.list <- lapply(file.list, "read.table", sep = '\t', stringsAsFactors = F, header = T)
names(df.list) <- str_replace_all(file.list, "(^out/)|(\\.txt$)", "")

# filtered df
# filter.dir <- "out/2019-07-30_directional_phenotype_scores_all.txt"
# filter.undir <- "out/2019-07-30_undir_scores_all.txt"
# filter.dir.df <- read.delim(unfilter.dir, sep = "\t", stringsAsFactors = F)
# filter.undir.df <- read.delim(unfilter.undir, sep = '\t', stringsAsFactors = F)
# 

# setup
eom <- paste0("GAMuT_", read.csv("EOMexclude.csv", header = F, stringsAsFactors = F)[,1])
plotScore <- function(df, value = "TotalScore", exl = eom, fn = NULL) {
  mat <- dcast(df, phenotype ~ GAMuT_ID, value.var = value)
  # rownames(mat) <- mat$phenotype
  rownames(mat) <- mat$phenotype
  mat$phenotype <- NULL
  mat <- as.matrix(mat)
  mat <- mat[order(rowSums(mat), decreasing = T), ]
  
  column.col <- ifelse(colnames(mat) %in% exl, "green", "red")
  if (is.null(fn)) {
    fn <- paste(Sys.Date(), deparse(substitute(df)), significance, "totalScore.pdf", sep = "_")
  }
  pdf(file = fn, onefile = T, width = 10, height = 10)

  heatmap.2(mat, ColSideColors = column.col, labRow = F, 
            dendrogram = "column", key = T, margins = c(7, 7), density.info = "none", 
            trace = "none", scale = "none", keysize = 1.5, key.ylab = "Total score")
  
  dev.off()
}

plotScore2 <- function(df, sig.count, labels, fn = NULL, 
                       thres = 0.01, probs = 0.5) {
  # takes into account significance
  # filter for: 1. no ben show significance, 2. 0.75% percentile cut-off for moc (top 0.25)
  sig.ben <- sig.count$ben.sig == 0
  sig.moc <- sig.count$moc.sig > quantile(sig.count$moc.sig, probs = probs)
  sig.phe <- sig.count[sig.ben & sig.moc, "phenotype"]
  df <- df[df$phenotype %in% sig.phe, ]
  df <- df[df$GAMuT_ID %in% rownames(labels)[labels$Labels == "MOC"], ]
  df[df$p.val >= 0.05, "TotalScore"] <- -100
  # df[, "is.sig"] <- df$p.val < 0.05
  mat <- dcast(df, phenotype ~ GAMuT_ID, value.var = "TotalScore")
  rownames(mat) <- mat[,1]
  mat <- mat[,-1]
  pheno.clst <- hclust(dist(mat))
  samp.clst <- hclust(dist(t(mat)))
  
  # build dendrograms
  pheno.dendr <- dendro_data(pheno.clst)
  samp.dendr <- dendro_data(samp.clst)
  
  pheno.dg <- ggplot(data = segment(pheno.dendr)) + 
    geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) + 
    coord_flip()
  
  samp.dg <- ggplot(data = segment(samp.dendr)) + 
    geom_segment(aes(x=x,y=y,xend=xend,yend=yend))
  
  # reorder the columns and rows
  df$phenotype <- factor(df$phenotype, levels = pheno.clst$labels[pheno.clst$order])
  df$GAMuT_ID <- factor(df$GAMuT_ID, levels = samp.clst$labels[samp.clst$order])
  df[df$p.val >= 0.05, "TotalScore"] <- NA
  
  hmp <- ggplot(data = df) + 
    geom_tile(aes(x = GAMuT_ID, y = phenotype, fill = TotalScore)) + 
    scale_fill_gradient(low = "red", high = "yellow", na.value = "grey60") + 
    theme(axis.text.x = element_text(angle=90))
  
  pheno.strw <- max(strwidth(df$phenotype, units = 'inches'))
  samp.strw <- max(strwidth(df$GAMuT_ID, units = 'inches'))
  
  return(list(phe = pheno.dg, samp = samp.dg, hmp = hmp, 
              pheno.strw = pheno.strw, samp.strw = samp.strw))
}

comb.plotscore <- function(plot.res, top = "samp", main = "hmp", right = "phe", 
                           tr = 0.07, tl = 3.9, rt = 0.07, rb = 1.37) {
  top <- plot.res[[top]]
  main <- plot.res[[main]]
  right <- plot.res[[right]]
  
  main.grob <- ggplot_gtable(ggplot_build(main))$grobs
  lgd <- main.grob[[which(sapply(main.grob, function(e) e$name == "guide-box"))]]
  
  dg.theme <- theme(axis.title = element_blank(), axis.ticks = element_blank(), 
                    axis.text = element_blank(), panel.background = element_blank(), 
                    panel.grid = element_blank()) 
  
  top <- top + dg.theme + 
    theme(plot.margin = unit(c(0, tr, 0, tl), 
                             unit = "inches")) + 
    scale_y_continuous(expand = expand_scale()) + 
    scale_x_continuous(expand = expand_scale(add = c(0.6, 0.6)))
    
  right <- right + dg.theme + 
    theme(plot.margin = unit(c(rt, 0, rb, 0), 
                             unit = "inches")) + 
    scale_y_continuous(expand = expand_scale()) + 
    scale_x_continuous(expand = expand_scale(add = c(0.6, 0.6)))
  
  main.clean <- main + theme(legend.position = "none", 
                             axis.text.y = element_text(size = 8)) + 
    scale_y_discrete(expand = expand_scale()) + 
    scale_x_discrete(expand = expand_scale()) 
  grid.arrange(top, lgd, main.clean, right, ncol = 2, nrow = 2, 
               heights = unit(c(3, 8+plot.res$samp.strw), units = "inches"), 
               widths = unit(c(6+plot.res$pheno.strw, 3), units = "inches"))

}

plotScore3 <- function(df, sig.count, labels, el, annotation, annot.col, fn = NULL, 
                       thres = 0.01, probs = 0.15, method = "euclidean", 
                       cluster_method = "complete", ...) {
  
  sig.ben <- sig.count$ben.sig == 0
  sig.moc <- sig.count$moc.sig > probs * length(which(labels$Labels == "MOC"))
  sig.phe <- sig.count[sig.ben & sig.moc, "phenotype"]
  df <- df[df$phenotype %in% sig.phe, ]
  df <- df[df$GAMuT_ID %in% rownames(labels)[labels$Labels == "MOC"], ]
  sig.row <- df$p.val < thres
  # create matrix for plotting, non sig values as NA
  
  for.plot <- dcast(df[sig.row, ], phenotype ~ GAMuT_ID, value.var = "TotalScore", fill = NA)
  rownames(for.plot) <- for.plot[, 1]
  for.plot <- as.matrix(for.plot[, -1])
  
  df[!sig.row, "TotalScore"] <- -100
  # df[!sig.row & df$TotalScore < 0, "TotalScore"] <- -100
  # df[!sig.row & df$TotalScore >= 0, "TotalScore"] <- 100
  # create matrix for clustering (non sig values infinitely far)
  keep1 <- df$GAMuT_ID %in% colnames(for.plot)
  keep2 <- df$phenotype %in% rownames(for.plot)
  for.clst <- dcast(df[keep1 & keep2, ], phenotype ~ GAMuT_ID, value.var = "TotalScore")
  rownames(for.clst) <- for.clst[, 1]
  for.clst <- as.matrix(for.clst[,-1])
  
  row.clst <- hclust(dist(for.clst, method = method), method = cluster_method)
  col.clst <- hclust(dist(t(for.clst), method = method), method = cluster_method)
  
  # get subset annotations
  annotation <- annotation[col.clst$labels[col.clst$order], ]
  for.plot <- for.plot[rownames(for.clst), colnames(for.clst)]
  
  p <- pheatmap(for.plot, cluster_rows = row.clst, cluster_cols = col.clst,
                annotation_col = annotation, annotation_colors = annot.col, ...)
  
  return(p)
  
}


plot.dir <- paste0(Sys.Date(), "_plot/")
if (!(plot.dir %in% dir())) dir.create(plot.dir)

# plot heatmaps of the score df's
for (n in names(df.list)) {
  plotScore(df.list[[n]], fn = paste0(plot.dir, n, "_plotScore.pdf"))
  plotScore2(df.list[[n]], fn = paste0(plot.dir, n, "_plotScore2.pdf"))
}

# plot dispersion and rank density for each significant phenotype
# all.sig <- getAllSignificantPhenotype(df.list, thres = 0.01)
all.sig <- getBenNegative(df.list, thres = 0.01, lb = labels)
for (n in names(all.sig)) {
  paired <- grepl("_paired", n)
  df <- df.list[[n]]
  for (p in all.sig[[n]]) {
    if (paired) {
      gs.up <- ptw[[paste0(p, "_UP")]]
      gs.dn <- ptw[[paste0(p, "_DN")]]
    } else {
      gs <- ptw[[p]]
    }
    if (!(p %in% dir(plot.dir))) {
      dir.create(paste0(plot.dir, p))
    }
    fn <- paste0(plot.dir, p, "/DispersionPlot.pdf")
    pdf(file = fn, width = 9, height = 7)
    # plot dispersion
    disp.df <- df[df$phenotype == p, setdiff(colnames(df), c("GAMuT_ID", "phenotype", "p.val", "Type"))]
    plot(plotDispersion(disp.df, annot = labels$Labels))
    dev.off()
    
    # plot rank density for each patient
    for (s in colnames(rank.data)) {
      fn <- paste0(plot.dir, p, "/", paste(s, "RankDensity.pdf", sep = "_"))
      pdf(file = fn, width = 9, height = 7)
      if (paired) {
        plot(plotRankDensity(rankData = rank.data[, s, drop = F], upSet = gs.up, 
                        downSet = gs.dn, isInteractive = F))
        
      } else {
        plot(plotRankDensity(rankData = rank.data[, s, drop = F], upSet = gs,
                        isInteractive = F))
      }
      dev.off()
    }
  }
}

ca125 <- paste0("GAMuT_", 
                read.csv("../CA125exclude.csv", header = F, stringsAsFactors = F)[,1])
df.non.ca <- lapply(df.list, function(e) {
  return(e[!(e$GAMuT_ID %in% ca125), ])
})

# example code: plotting landscape -----
df1 <- df.list$`2019-08-05_paired.dir_c2.all_FALSE_results`
df2 <- df.list$`2019-08-05_unpaired.dir_c2.all_FALSE_results`
df.liver <- df1[df1$phenotype == "KIM_LIVER_CANCER_POOR_SURVIVAL",]
df.liver2 <- df2[df2$phenotype == "IIZUKA_LIVER_CANCER_PROGRESSION_L0_L1_DN",]
rownames(df.liver) <- df.liver$GAMuT_ID
rownames(df.liver2) <- df.liver2$GAMuT_ID

# plotScoreLandscape is a plot of one gene set vs. another for comparison
# projectScoreLandscape annotate the plot produced by the previous function 
liver.plot <- plotScoreLandscape(df.liver, df.liver2, 
                                scorenames = c('liver','liver2'), isInteractive = FALSE)

projectScoreLandscape(plotObj = liver.plot, df.liver, df.liver2, 
                      subSamples = rownames(df.liver)[1:4], 
                      annot = rownames(df.liver)[1:4], 
                      sampleLabels = round(df.liver$TotalScore[1:4],2),isInteractive = F)

# plotDispersion plots the score(total, up, and down) against dispersion
# with colour codes to indicate different phenotype
disp.df <- df.liver[, setdiff(colnames(df.liver), c("phenotype", "GAMuT_ID", "p.val", "Type"))]
plotDispersion(disp.df, labels$Labels, isInteractive = F)

# plotRankDensity plots the running sum plots similar to the GSEA plot
# this for each patient 
liver.up <- ptw[["KIM_LIVER_CANCER_POOR_SURVIVAL_UP"]]
liver.dn <- ptw[["KIM_LIVER_CANCER_POOR_SURVIVAL_DN"]]
plotRankDensity(rankData = rank.data[, 20, drop = F], upSet = liver.up, downSet = liver.dn, isInteractive = F)

colnames(rank.data) <- str_replace(colnames(rank.data), "^GAMuT_", "")
rownames(labels) <- str_replace(rownames(labels), "^GAMuT_", "")

tmp <- df.list$`2019-08-05_paired.dir_c2.all_FALSE_results`
tmp2 <- sig.count.list$`out/2019-08-05_paired.dir_c2.all_FALSE_results.txt`
keep.phe <- tmp2[tmp2$ben.sig == 0, "phenotype"]
tmp <- tmp[tmp$p.val < 0.05 & tmp$phenotype %in% keep.phe & 
             labels[tmp$GAMuT_ID, ] == "MOC", ]
tmp <- tmp[order(tmp$TotalScore, decreasing = T), ]
for (i in 1:6) {
  # paired phenotypes
  r <- tmp[i, ]
  up.set = ptw[[paste0(r$phenotype, "_UP")]]
  dn.set = ptw[[paste0(r$phenotype, "_DN")]]
  fn <- paste(r$phenotype, r$GAMuT_ID, r$TotalScore, "pdf", sep = ".")
  pdf(height = 5, width = 7, file = paste0("2019-10-02_plot/paired-dir/", fn))
  grid.arrange(plotRankDensity(rankData = rank.data[, r$GAMuT_ID, drop = F], 
                                       upSet = up.set, downSet = dn.set,
                                        isInteractive = F))
  dev.off()
}
# unpaired directional
tmp <- df.list$`2019-08-05_unpaired.dir_c2.all_FALSE_results`
tmp2 <- sig.count.list$`out/2019-08-05_unpaired.dir_c2.all_FALSE_results.txt`
keep.phe <- tmp2[tmp2$ben.sig == 0, "phenotype"]
tmp <- tmp[tmp$p.val < 0.05 & tmp$phenotype %in% keep.phe & 
             labels[tmp$GAMuT_ID, ] == "MOC", ]
tmp <- tmp[order(tmp$TotalScore, decreasing = T), ]

for (i in 1:6) {
  r <- tmp[i, ]
  up.set <- ptw[[r$phenotype]]
  fn <- paste(r$phenotype, r$GAMuT_ID, r$TotalScore, "pdf", sep = ".")
  pdf(height = 5, width = 7, file = paste0("2019-10-02_plot/unpaired-dir/", fn))
  grid.arrange(plotRankDensity(rankData = rank.data[, r$GAMuT_ID, drop = F], 
                               upSet = up.set, isInteractive = F))
  dev.off()
}

# non-directional
tmp <- df.list$`2019-08-05_undir_c2.all_FALSE_results`
tmp2 <- sig.count.list$`out/2019-08-05_undir_c2.all_FALSE_results.txt`
keep.phe <- tmp2[tmp2$ben.sig == 0, "phenotype"]
tmp <- tmp[tmp$p.val < 0.05 & tmp$phenotype %in% keep.phe & 
             labels[tmp$GAMuT_ID, ] == "MOC", ]
tmp <- tmp[order(tmp$TotalScore, decreasing = T), ]
for (i in 1:6) {
  r <- tmp[i, ]
  up.set <- ptw[[r$phenotype]]
  fn <- paste(r$phenotype, r$GAMuT_ID, r$TotalScore, "pdf", sep = ".")
  pdf(height = 5, width = 7, file = paste0("2019-10-02_plot/undir/", fn))
  grid.arrange(plotRankDensity(rankData = rank.data[, r$GAMuT_ID, drop = F], 
                               upSet = up.set, isInteractive = F))
  dev.off()
}


## read result files and determine which pathway is significant not in ben but in MOC
file <- paste0("out/", dir(path = "out/", pattern = "^2019-08-05.*results.txt"))
moc.ids <- rownames(labels)[labels$Labels == "MOC"]
ben.ids <- rownames(labels)[labels$Labels == "BEN"]
sig.count.list <- list()
for (f in file) {
  res.df <- read.delim(f, header = T, sep = "\t", stringsAsFactors = F)
  chs.moc <- res.df[, "GAMuT_ID"] %in% moc.ids
  chs.ben <- res.df[, "GAMuT_ID"] %in% ben.ids
  pheno.n <- unique(res.df[, "phenotype"])
  sig.count <- lapply(pheno.n, FUN = function(n) {
    chs.phe <- res.df[, "phenotype"] == n
    moc.sig <- length(which(res.df[chs.moc & chs.phe, "p.val"] < 0.05))
    ben.sig <- length(which(res.df[chs.ben & chs.phe, "p.val"] < 0.05))
    return(data.frame(phenotype = n, moc.sig = moc.sig, ben.sig = ben.sig))
  })
  sig.count <- do.call("rbind", sig.count)
  sig.count[,1] <- as.character(sig.count[,1])
  sig.count.list[[f]] <- sig.count
}

for (sig.count in sig.count.list) {
  x <- sig.count[sig.count$ben.sig == 0, ]
}

# load cnv and variant data 
load("../GoldDataUseWhenEverPossible/mut.and.exp.matrix.RData")
getDriverEvent <- function() {}
insertSource("../GitHubDriverNet/PlottingFunctions.R", functions = "getDriverEvent")
getDriverEvent <- getDriverEvent@.Data

genes <- c("KRAS", "ERBB2", "TP53", "CDKN2A")
mocCNVData$Sample <- str_replace(mocCNVData$Sample, "^GAMuT_", "")
mocVariantData$Sample <- str_replace(mocVariantData$Sample, "^GAMuT_", "")

el <- getDriverEvent(genes, id.name = "Sample")
# cons.rank.labels = c("1" = "truncating frame-shift", 
#            "2" = "essential splice site", 
#            "3" = "in-frame indel", 
#            "4" = "stop codon loss", 
#            "NONE" = "None")


getAnnotation <- function(samples, cnv, var, el) {

  mat <- data.frame(row.names = samples)
  g.list = list(cnv = cnv, variant = var)
  for (kind in names(g.list)) {
    for (g in g.list[[kind]]) {
      cn <- paste(g, kind, sep = ".")
      val <- ifelse(kind == 'cnv', "Event", "Consequence_Rank")
      mat[, cn] <- sapply(samples, FUN = function(e) {
        tmp <- el[[g]][[kind]]
        v <- as.character(tmp[tmp$GAMuT_ID == e, val, drop = T])
        return(ifelse(length(v) == 0, "None", v))
      })
      
    }
  }
  
  for (i in grep("variant", colnames(mat))) {
    mat[, i] <- consequence.rank.conv(mat[, i])
  }
  
  for (i in grep("cnv", colnames(mat))) {
    mat[, i] <- factor(mat[, i], levels = c("None", "HA", "HD"))
  }
  return(mat)
}

consequence.rank.conv <- function(x) {
  guide <- c("1" = "truncating frame-shift", 
             "2" = "essential splice site",
             "3" = "in-frame indel",
             "4" = "stop codon loss")
  lvs <- c("None", guide)
  x <- as.character(x)
  for (rank in names(guide)) {
    x[x == rank] <- guide[rank]
  }
  x <- factor(x, levels = lvs)
  return(x)
}

getAnnotColour <- function(annot) {
  annot.colour <- list()
  annot.colour[["Stage"]] <- c("I" = "palegreen", "II" = "yellow", "III-IV" = "red")
  annot.colour[["CA125.high"]] <- c("No" = "grey90", "Yes" = "grey40")
  annot.colour[["Grade"]] <- c("1" = "steelblue1", "2" = "steelblue3", "3" = "steelblue4")
  i <- grep("cnv$", colnames(annot), value = T)
  annot.colour[i] <- list(c("None" = "grey90", "HD" = "steelblue2", "HA" = "firebrick1"))
  i <- grep("variant$", colnames(annot), value = T)
  annot.colour[i] <- list(c("None" = "grey90", "truncating frame-shift" = "red", 
                       "essential splice site" = "pink", "in-frame indel" = "purple",
                       "stop codon loss" = "orange"))
  return(annot.colour)
}


# remove the GAMuT_ prefix
for (i in 1:length(df.list)) {
  df.list[[i]]$GAMuT_ID <- str_replace(df.list[[i]]$GAMuT_ID, "^GAMuT_", "")
}
rownames(labels) <- str_replace(rownames(labels), "^GAMuT_", "")

annot.row <- getAnnotation(samples = unique(df.list[[1]]$GAMuT_ID), 
                          cnv = c("ERBB2"), var = c("KRAS", "TP53", "CDKN2A"), el = el)

ca125 <- str_replace(ca125, "^GAMuT_", "")

stage.data <- read.csv("All survival_CN_Aug19.csv", as.is = T)
rownames(stage.data) <- stage.data$GAMUT_ID
annot.row$Stage <- stage.data[rownames(annot.row), "Stage.group.cat"]
annot.row["844.07G3", "Stage"] <- "I"
annot.row["OV2602.A", "Stage"] <- "III-IV"

annot.row$Grade <- stage.data[rownames(annot.row), "Grade"]
annot.row["844.07G3", "Grade"] <- "3"
annot.row["OV2602.A", "Grade"] <- "2"
annot.row$Grade <- factor(annot.row$Grade, levels = c(1, 2, 3))


annot.row$CA125.high <- "No"
annot.row[rownames(annot.row) %in% ca125, "CA125.high"] <- "Yes"
annot.row$CA125.high <- factor(annot.row$CA125.high, levels = c("No", 'Yes'))

##### found EOM that was found to be MOC -------
new.eom <- c("22826")
for (i in 1:length(df.list)) {
  df.list[[i]] <- df.list[[i]][!(df.list[[i]]$GAMuT_ID %in% new.eom), ]
}
labels <- labels[!(rownames(labels) %in% new.eom), , drop = F]
stage.lvs <- c("I", "II", "III-IV")
annot.row <- annot.row[annot.row$Stage %in% stage.lvs, ]
annot.row$Stage <- factor(annot.row$Stage, levels = stage.lvs)
annot.colour <- getAnnotColour(annot.row)

color <- colorRampPalette(c("red", "yellow"))(100)

plotScore3(df = df.list[[1]], sig.count = sig.count.list[[1]], labels = labels, el = el, 
           probs = 0.15, 
           annotation = annot.row, annot.col = annot.colour, thres = 0.05, method = "manhattan", 
           cluster_method = "ward.D2", color = color, cellwidth = 12, cellheight = 12, 
           treeheight_row = 120, treeheight_col = 120, 
           filename = "2019-10-11_paired.dir_hmp_pval_0.05.pdf")

plotScore3(df = df.list[[2]], sig.count = sig.count.list[[2]], labels = labels, el = el, 
           probs = 0.15, 
           annotation = annot.row, annot.col = annot.colour, thres = 0.05, method = "manhattan", 
           cluster_method = "ward.D2", color = color, cellwidth = 12, cellheight = 12, 
           treeheight_row = 120, treeheight_col = 120, 
           filename = "2019-10-11_undir_hmp_pval_0.05.pdf")

plotScore3(df = df.list[[3]], sig.count = sig.count.list[[3]], labels = labels, el = el, 
           probs = 0.15, 
           annotation = annot.row, annot.col = annot.colour, thres = 0.05, method = "manhattan", 
           cluster_method = "ward.D2", color = color, cellwidth = 12, cellheight = 12, 
           treeheight_row = 120, treeheight_col = 120, 
           filename = "2019-10-11_unpaired.dir_hmp_pval_0.05.pdf")

dev.off()

# get significant phenotypes in each group
sig.phe <- list()
for (n in names(sig.count.list)) {
  sig.phe[[n]] <- sig.count.list[[n]][sig.count.list[[n]]$ben.sig == 0 & sig.count.list[[n]]$moc.sig > 0, 
                                      "phenotype"]
}
# load GSEA results for comparison
load("../GSEA/2019-10-11_sig.ptw.RData")

names(sig.phe)
names(sig.phe) <- c("paired.dir", "undir", "unpaired.dir")
sig.phe[["paired.dir.sep"]] <- c(paste0(sig.phe[["paired.dir"]], "_UP"), 
                                 paste0(sig.phe[["paired.dir"]], "_DN"))

save(list = "sig.phe", file = "2019-10-13_sig.phe.RData")
#see intersections
comm.ptw <- list()
for (n in names(sig.phe)[2:4]) {
  comm.ptw[[n]] <- intersect(sig.phe[[n]], names(sig.ptw))
}

length(unlist(comm.ptw))
# [1] 27
# total 27 ptws in all groups

save(list = "comm.ptw", file = "2019-10-13_comm.ptw.RData")

# test association of scores with stage, ca125, and grade
# a lot of directional pathwas also have a "non-directional" version
# e.g. "KYNG_RESPONSE_TO_H2O2_VIA_ERCC6_UP" "KYNG_RESPONSE_TO_H2O2_VIA_ERCC6_DN"
# "KYNG_RESPONSE_TO_H2O2_VIA_ERCC6"
# will add a suffix to differentiate them

df.list.cat <- lapply(df.list, FUN = function(e) {
  return(e[, c('GAMuT_ID', 'TotalScore', 'TotalDispersion', 'phenotype', 'p.val', 'Type')])
})
names(df.list.cat) 
names(df.list.cat) <- c("paired.dir", "undir", "unpaired.dir")
# add suffix to differentiate similar named pathways
# also filter for significant pathways
moc.ids <- str_replace(moc.ids, "^GAMuT_", "")
for (n in names(df.list.cat)) {
  df.list.cat[[n]] <- df.list.cat[[n]][df.list.cat[[n]]$phenotype %in% sig.phe[[n]] &
                                         df.list.cat[[n]]$GAMuT_ID %in% moc.ids, ]
  df.list.cat[[n]]$phenotype <- paste0(df.list.cat[[n]]$phenotype, "_", n)
}
df.cat <- do.call("rbind", df.list.cat)
# fill by NA neglects non-significant samples, so only significant scores were kept
df.cat.dmat <- dcast(df.cat, GAMuT_ID ~ phenotype, value.var = "TotalScore", fill = NA)
dim(df.cat.dmat)
rownames(df.cat.dmat) <- df.cat.dmat$GAMuT_ID
test.phe <- colnames(df.cat.dmat)[-1]
df.cat.dmat <- cbind(df.cat.dmat, annot.row[rownames(df.cat.dmat), ])
# only two stage II patients, combining them with stage III-IV
df.cat.dmat$Stage <- as.character(df.cat.dmat$Stage)
df.cat.dmat$Stage[df.cat.dmat$Stage == "II"] <- "II-IV"
df.cat.dmat$Stage[df.cat.dmat$Stage == "III-IV"] <- "II-IV"
df.cat.dmat$Stage <- factor(df.cat.dmat$Stage, levels = c("I", "II-IV"))

pval.mat <- matrix(NA, length(test.phe), 3, dimnames = list(test.phe, c("Stage", "CA125.high", "Grade")))
for (fact in c("Stage", "CA125.high", "Grade")) {
  for (p in test.phe) {
    fml.str <- paste(p, " ~ ", fact)
    fml <- as.formula(fml.str)
    aov.res <- aov(fml, data = df.cat.dmat)
    pval.mat[p, fact] <- unlist(summary.aov(aov.res))["Pr(>F)1"]
  }
}
write.table(as.data.frame(pval.mat),file = "2019-10-13_anova_singscore.txt", sep = "\t",
            col.names = T, row.names = T, quote = F)

pval.melt <- melt(pval.mat, varnames = c("Phenotype", "Group"), value.name = "p.value")
pval.melt$Phenotype <- as.character(pval.melt$Phenotype)
pval.melt$Group <- as.character(pval.melt$Group)
pval.melt <- pval.melt[order(pval.melt$p.value), ]

# plot for top n most significant pairs for each group
n <- 6
top.tab <- NULL
for (fact in c("Stage", "CA125.high", "Grade")) {
  tmp.df <- pval.melt[pval.melt$Group == fact, ][1:n, ]
  top.tab <- rbind(top.tab, tmp.df)
  for (i in 1:n) {
    p <- tmp.df[i, "Phenotype"]
    plot.title <- str_replace(p, "_(paired.dir|undir|unpaired.dir)$", "")
    dt <- df.cat.dmat[, c(p, fact)]
    pdf(file = paste0("2019-10-13_dotplot_anova/", paste(p, fact, "pdf", sep = ".")), height = 3, width = 4)
    bw <- (max(dt[, p]) - min(dt[, p])) / 50
    grid.arrange(ggdotplot(data = dt, x = fact, y = p, color = fact, ylab = "Score", xlab = fact,
                           add = c("mean_se", 'boxplot'), binwidth = bw, error.plot = "errorbar",
                           size = 2, title = plot.title) +
                   theme(legend.position = "none", text = element_text(size = 6),
                         plot.title = element_text(size = 8, hjust = 0)))
    dev.off()
  }
}
write.table(top.tab, file = "2019-10-13_top_significant_association.txt", sep = "\t", 
            quote = F, col.names = T, row.names = F)


comm.phe <- comm.ptw
comm.phe$paired.dir <- unique(str_replace(comm.phe$paired.dir.sep, "_(DN|UP)$", ""))
comm.phe <- comm.phe[names(comm.phe) != "paired.dir.sep"]

sig.count.all <- do.call(rbind, sig.count.list)
moc.sig.all <- data.frame(count = sig.count.all[sig.count.all$phenotype %in% unlist(sig.phe), "moc.sig"],
                          group = "All")
moc.sig.comm <- data.frame(count = sig.count.all[sig.count.all$phenotype %in% unlist(comm.phe), "moc.sig"], 
                           group = "Common")
moc.sig.dt <- rbind(moc.sig.all, moc.sig.comm)
ggboxplot(data = moc.sig.dt, x = "group", y = "count", colour = "count", 
          add = "boxplot")

gghistogram(sig.count.all[sig.count.all$phenotype %in% unlist(comm.phe), "moc.sig"], 
            bins = max(sig.count.all$moc.sig), fill = "black")
gghistogram(sig.count.all[sig.count.all$phenotype %in% unlist(sig.phe), "moc.sig"], 
            bins = max(sig.count.all$moc.sig), fill = "black")
