rm(list =ls())
setwd("/Volumes/GoogleDrive/My Drive/DWV_TRN/HB_NewerAnnotation/Count_tables/prefiltered.universe/")
library(edgeR)
library(limma)

counts = read.delim("Hb_challenge_drone_ftzf1_gbb.2019-01-12", header = T, sep = " ")
DGE_all = DGEList(counts=counts)

#Expression = cpm(DGE_all, log = T)
Expression = cpm(DGE_all, log = T)

setwd("/Volumes/GoogleDrive/My Drive/DWV_TRN/")
hb_tfs = read.delim("hb_tfs.txt", header = F, sep = "\t")
hb_tfs_expressed = intersect(row.names(Expression), hb_tfs$V1)

write.table(hb_tfs_expressed, file = paste("Hb_scDWV.eTFs.txt", Sys.Date(), "txt",sep = "."), quote = F, col.names = F, row.names = F)
write.table(Expression, file = paste("Hb_scDWV.Exprs.Astrix", Sys.Date(), "txt",sep = "."), sep = "\t", quote = F, col.names = F, row.names = F)
write.table(rownames(Expression), file = paste("Hb_scDWV.Probes.Astrix", Sys.Date(), "txt",sep = "."), sep = "\t", quote = F, col.names = F, row.names = F)
