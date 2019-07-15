rm(list =ls())
setwd("/Users/syedabbasbukhari/Google Drive/social_challenge/honey_DMV_TRN/")
library(edgeR)
library(limma)
library(preprocessCore)
targets = readTargets("hb_targets.txt")


getDGE = function(targets)
{
  Raw_DGE = readDGE(targets)
  MetaTags <- grep("^__", rownames(Raw_DGE))
  
  Raw_DGE = Raw_DGE[-MetaTags,]
  
  keep <- rowSums(cpm(Raw_DGE)>0.5) >= 5
  Filtered_DGE = Raw_DGE[keep,]
  Filtered_DGE = calcNormFactors(Filtered_DGE)
}

### Getting TRN ###

DGE_all = getDGE(targets)

#Expression = cpm(DGE_all, log = T)
Expression = DGE_all$counts
Expression = normalize.quantiles(Expression)
Expression <- scale(Expression)
hb_tfs = read.delim("hb_tfs.txt", header = F, sep = "\t")
hb_tfs_expressed = intersect(row.names(Expression), hb_tfs$V1)

write.table(hb_tfs_expressed, file = paste("Hb_scDWV.eTFs.txt", Sys.Date(), "txt",sep = "."), quote = F, col.names = F, row.names = F)
write.table(Expression, file = paste("Hb_scDWV.Exprs.Astrix", Sys.Date(), "txt",sep = "."), sep = "\t", quote = F, col.names = F, row.names = F)
write.table(rownames(Expression), file = paste("Hb_scDWV.Probes.Astrix", Sys.Date(), "txt",sep = "."), sep = "\t", quote = F, col.names = F, row.names = F)
