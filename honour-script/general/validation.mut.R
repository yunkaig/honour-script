# use expanded mutation data (patients without RNA-seq))

all.var <- read.csv("mutExpData/allVariantData.csv", stringsAsFactors = F)
old.cnv <- read.delim("mutExpData/allCNVData.txt", stringsAsFactors = F)
old.cnv <- old.cnv[old.cnv$Event %in% c("High Copy Gain", "Homozygous Copy Loss"), ]
all.var <- all.var[all.var$Consequence_Rank < 5, ]

length(unique(old.cnv$Sample))

# exclude Sanger Key genes Samples (See supplementary table)
# split by source
# PCA plot of RNA-seq data (see sensitivity)


## load Name.csv file for variant sample ids
name.comp <- read.csv("mutExpData/Names.csv", stringsAsFactors = F)
name.comp2 <- read.csv("mutExpData/Name comparison.csv", stringsAsFactors = F)

# variant samples have GAMuT_ID within the file 

cnv.gamut <- vapply(old.cnv$Sample, FUN.VALUE = "1", FUN = function(e) {
  return(name.comp[which(name.comp$Nexus.ID == e), "GAMUT_ID"])
})

any(is.na(cnv.gamut))
old.cnv$GAMuT_ID <- cnv.gamut

write.table(old.cnv, "mutExpData/allCNVDataFilteredIDConverted.txt", quote = F, 
            col.names = T, row.names = F, sep = "\t")


## adding the new CNV data (April 2019 WM samples)
new.cnv <- read.delim("mutExpData/20190724_Nexus Export WM new samples.txt", stringsAsFactors = F)
new.cnv <- new.cnv[new.cnv$Event %in% c("High Copy Gain", "Homozygous Copy Loss"), ]
new.cnv.gamut <- vapply(new.cnv$Sample, FUN.VALUE = "", FUN = function(e) {
  if (!(e %in% name.comp$Nexus.ID)) return(as.character(NA))
  return(name.comp[which(name.comp$Nexus.ID == e), "GAMUT_ID"])
})

unique(new.cnv[is.na(new.cnv.gamut), "Sample"])
new.cnv.gamut[is.na(new.cnv.gamut)] <- name.comp[grep("WM432", name.comp$Nexus.ID), "GAMUT_ID"]
any(is.na(new.cnv.gamut))

new.cnv$GAMuT_ID <- new.cnv.gamut
write.table(new.cnv, "mutExpData/20190724_Nexus Export WM new samples Filtered ID converted.txt", quote = F, 
            col.names = T, row.names = F, sep = "\t")

## combine the new with the old cnv files
tmp <- str_split(new.cnv$Chromosome.Region, ":")
new.chr <- sapply(tmp, FUN = function(e) e[1])
new.region <- sapply(tmp, FUN = function(e) e[2])
new.cnv$Chromosome <- new.chr
new.cnv$Region <- new.region
new.cnv$Chromosome.Region <- NULL

cols.keep <- c("Sample", "GAMuT_ID", "Chromosome", "Region", "Event", "Length", "Cytoband", "Gene.Symbols")
all.cnv <- rbind(old.cnv[, cols.keep],
                  new.cnv[, cols.keep])

write.table(all.cnv, "mutExpData/combinedCNV_filtered_GAMUT.txt", quote = F, 
            col.names = T, row.names = F, sep = "\t")



# types of analysis run on samples; get valid CNVs and Variant samples
an.table <- read.csv("mutExpData/All survival_CN_Aug19.csv", stringsAsFactors = F)
fil1 <- an.table$Classification %in% c("MOC", "BDL", "BEN")      # filter for MOT cases
fil2 <- !(an.table$Sequencing.panel %in% c("Sanger key genes", "Failed"))    # filter for valid sequencing.panel
valid.samples <- an.table[fil1 & fil2, ]

# couldn't find the following samples in an.table$GAMUT_ID
#
# > setdiff(all.var$GAMuT_ID, an.table$GAMUT_ID)
# [1] "C1711"   "OV2602C" "WM982A" 

mocGamutIDs <- valid.samples[valid.samples$Classification == "MOC", "GAMUT_ID"]
bdlGamutIDs <- valid.samples[valid.samples$Classification == "BDL", "GAMUT_ID"]
benGamutIDs <- valid.samples[valid.samples$Classification == "BEN", "GAMUT_ID"]

length(mocGamutIDs)
# 243
length(bdlGamutIDs)
# 108
length(benGamutIDs)
# 5


valid.var <- all.var[all.var$GAMuT_ID %in% mocGamutIDs, ]
valid.cnv <- all.cnv[all.cnv$GAMuT_ID %in% mocGamutIDs, ]
length(unique(valid.var$GAMuT_ID))
length(unique(valid.cnv$GAMuT_ID))
samp.var <- unique(valid.var$GAMuT_ID)
samp.cnv <- unique(valid.cnv$GAMuT_ID)

length(setdiff(samp.var, samp.cnv))
# 10
length(setdiff(samp.cnv, samp.var))
# 0
# the difference between the cnv samples and the variant samples is large
length(intersect(samp.cnv, samp.var))
# 180

el.valid <- getDriverEvent(uni.genes, cnv = valid.cnv, var = valid.var)
plotEventCountBar(genes.plot, el = el.valid)
plotGenomeAbberationHeatmap(genes = genes.plot, el = el.valid, show.cnv = T)
uni.col <- list(tomato = intersect(dawn.genes, genes.plot), 
                palegreen2 = setdiff(genes.plot, dawn.genes),
                steelblue2 = setdiff(dawn.genes, genes.plot))
plotPDF(plotHeatmapAndBar(uni.genes, el = el.valid, group.col = uni.col, gene.lab = T),
        file = "2019-09-05_heatmap_dendr_unigenes.pdf")
length(setdiff(all.cnv$GAMuT_ID, all.var$GAMuT_ID))


# save(list = c("mocTumourExpression", "mocCNVData", 
#               "mocVariantData", "valid.cnv", "valid.var", 
#               "benExpression"), file = "../GoldDataUseWhenEverPossible/all.variant.and.cnv.RData")


