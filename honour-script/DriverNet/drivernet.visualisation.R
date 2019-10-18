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

el <- getDriverEvent(genes.plot, cnv = valid.cnv, var = valid.var, id.name = "GAMuT_ID")
event.df <- getEventDF(uni.el)

# remove GAMUT prefix 
event.df$GAMuT_ID <- str_replace(event.df$GAMuT_ID, "^GAMuT_", "")


# get annotation from stage data
stage.data <- read.csv("../All survival_CN_Aug19.csv", as.is = T)
rownames(stage.data) <- stage.data$GAMUT_ID
annot.row <- data.frame(row.names = unique(event.df$GAMuT_ID))
annot.row$Stage <- stage.data[rownames(annot.row), "Stage.group.cat"]
annot.row["844.07G3", "Stage"] <- "I"
annot.row["OV2602.A", "Stage"] <- "III-IV"

# get annotation for grades
annot.row$Grade <- stage.data[rownames(annot.row), "Grade"]
annot.row["844.07G3", "Grade"] <- "3"
annot.row["OV2602.A", "Grade"] <- "2"
annot.row$Grade <- factor(annot.row$Grade, levels = c(1, 2, 3))

# remove new EOM case 
eom <- rownames(annot.row)[annot.row$Stage == "EOM"]
annot.row <- annot.row[setdiff(rownames(annot.row), eom), , drop = F]
event.df <- event.df[!(event.df$GAMuT_ID %in% eom), ]

ca125 <- read.csv("../CA125exclude.csv", as.is = T)[,1]
annot.row$CA125.high <- ifelse(rownames(annot.row) %in% ca125, "Yes", "No")

annot.colour <- list(Stage = c("I" = "palegreen", "II" = "yellow", "III-IV" = "tomato"), 
                     CA125.high = c("No" = "grey80", "Yes" = "grey40"),
                     Grade = c("1" = "steelblue1", "2" = "steelblue3", "3" = "steelblue4"))
plotMutation(df = event.df, annotation = annot.row, annot.col = annot.colour, 
             breaks = c(0.5, 1.5, 2.5, 3.5), color = c("steelblue2", "orangered1", "purple"), 
             cellheight = 12, cellwidth = 12, filename = "2019-10-11_hmp_unigenes_stage_ca125_grade_mut.pdf",
             cluster_method = "ward.D2")

event.df$Event.char <- as.character(event.df$Event)
all.pat <- rownames(annot.row)

event.dmat <- dcast(event.df, GAMuT_ID ~ SYMBOL, value.var = "Event.char", fill = "none")
rownames(event.dmat) <- event.dmat$GAMuT_ID
event.dmat <- event.dmat[rownames(annot.row), ]
event.dmat <- cbind(event.dmat, annot.row)

pval.mat <- matrix(NA, 3, length(uni.genes), 
                   dimnames = list(c("Stage", "CA125.high", "Grade"), uni.genes))
for (g in uni.genes) {
  stage.tab <- table(event.dmat$Stage, event.dmat[, g])
  pval.mat["Stage", g] <- fisher.test(stage.tab)$p.value
  ca125.tab <- table(event.dmat$CA125.high, event.dmat[, g])
  pval.mat["CA125.high", g] <- fisher.test(ca125.tab)$p.value
  grade.tab <- table(event.dmat$Grade, event.dmat[, g])
  pval.mat["Grade", g] <- fisher.test(grade.tab)$p.value
  
}
write.table(as.data.frame(pval.mat),file = "2019-10-11_fishers.test.uni.genes.txt", sep = "\t",
            col.names = T, row.names = T, quote = F)

# validate p-values with larger cohort, beware of the cohort with no RNA seq data

stage.data <- read.csv("../All survival_CN_Aug19.csv", as.is = T)
rownames(stage.data) <- stage.data$GAMUT_ID
annot.row <- data.frame(row.names = unique(valid.df$GAMuT_ID))
annot.row$Stage <- stage.data[rownames(annot.row), "Stage.group.cat"]
annot.row["844.07G3", "Stage"] <- "I"
annot.row["OV2602.A", "Stage"] <- "III-IV"
annot.row$Stage <- factor(annot.row$Stage, levels = c("I", "II", 'III-IV'))

# get annotation for grades
annot.row$Grade <- stage.data[rownames(annot.row), "Grade"]
annot.row["844.07G3", "Grade"] <- "3"
annot.row["OV2602.A", "Grade"] <- "2"
annot.row$Grade <- factor(annot.row$Grade, levels = c(1, 2, 3))

annot.row$CA125.high <- ifelse(rownames(annot.row) %in% ca125, "Yes", "No")
annot.row$CA125.high <- factor(annot.row$CA125.high, levels = c("No", 'Yes'))

valid.el <- getDriverEvent(genes = uni.genes, cnv = valid.cnv, var = valid.var)
valid.df <- getEventDF(valid.el)
valid.df$Event.char <- as.character(valid.df$Event)

plotMutation(df = valid.df, annotation = annot.row, annot.col = annot.colour, 
             breaks = c(0.5, 1.5, 2.5, 3.5), color = c("steelblue2", "orangered1", "purple"), 
             cellheight = 12, cellwidth = 12, 
             filename = "2019-10-11_hmp_unigenes_stage_ca125_grade_mut_validate.pdf",
             cluster_method = "ward.D2")

valid.dmat <- dcast(valid.df, GAMuT_ID ~ SYMBOL, value.var = "Event.char", fill = "none")

rownames(valid.dmat) <- valid.dmat$GAMuT_ID
valid.all.pat <- rownames(annot.row)
non.mut <- setdiff(valid.all.pat, rownames(valid.dmat))
# fix logic 2019-10-12: should use all patients with complete data
# and NOT ignore patients with no mutations in any driver genes

if (length(non.mut) > 0) {
  valid.dmat[non.mut, ] <- "none"
}
annot.row <- annot.row[rownames(valid.dmat), ]
valid.dmat <- cbind(valid.dmat, annot.row)

valid.rna.cohort <- intersect(rna.seq.cohort, rownames(valid.dmat))
pval.valid.mat <- matrix(NA, length(uni.genes), 3,
                         dimnames = list(uni.genes, c("Stage", "CA125.high", "Grade")))
for (g in uni.genes) {
  stage.tab <- table(valid.dmat[, "Stage"], valid.dmat[, g])
  pval.valid.mat[g, "Stage"] <- fisher.test(stage.tab)$p.value
  ca125.tab <- table(valid.dmat[valid.rna.cohort, "CA125.high"], 
                     valid.dmat[valid.rna.cohort, g])
  pval.valid.mat[g, "CA125.high"] <- fisher.test(ca125.tab)$p.value
  grade.tab <- table(valid.dmat[, "Grade"], valid.dmat[, g])
  pval.valid.mat[g, "Grade"] <- fisher.test(grade.tab)$p.value
}

write.table(as.data.frame(pval.valid.mat),file = "2019-10-12_fishers.test.uni.genes.validate.txt", 
            sep = "\t", col.names = T, row.names = T, quote = F)

#visualise top.genes
top.el <- getDriverEvent(top.genes, cnv = valid.cnv, var = valid.var, id.name = "GAMuT_ID")
plotEventCountBar(genes = top.genes, el = top.el, gene.lab = F)
plotHeatmapAndBar(top.genes, el = top.el)
