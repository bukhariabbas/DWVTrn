rm(list = ls())
setwd("/Volumes/GoogleDrive/My Drive/DWV_TRN/HB_NewerAnnotation/Count_tables/prefiltered.universe")

stopifnot(require(plyr))
stopifnot(require(readxl))
challenge =  read.delim("challenge.prefiltered.universe", header = T)
challenge$gid = rownames(challenge)
drone = read.delim("dronePER.prefiltered.universe", header = T)
drone$gid = rownames(drone)
ftzf1 = read.delim("ftzf1.prefiltered.universe", header = T)
ftzf1$gid = rownames(ftzf1)
gbb = read.delim("gbb.dance.prefiltered.universe", header = T)
gbb$gid = rownames(gbb)

combined = join(challenge, drone)
combined = join(combined, ftzf1)
combined = join(combined, gbb)
#adding gid as rownames
rownames(combined) = combined$gid
#dropping gid
combined = combined[, !(names(combined) %in% "gid")]

Annotation = as.data.frame(read_xlsx("HAv3.1_conversion_table.xlsx"))

Geneid_GBid = unique(Annotation$ID[which(Annotation$BEEBASE != "NA")])
combined_annotated = combined[which(rownames(combined) %in% Geneid_GBid),]

combined_annotated$gid = rownames(combined_annotated)
Geneid_GBid_pairs = unique(Annotation[which(Annotation$ID != "NA" & Annotation$BEEBASE != "NA"), c("ID", "BEEBASE")])
colnames(Geneid_GBid_pairs)[1] = "gid"
combined_annotated = join(combined_annotated, Geneid_GBid_pairs)
rownames(combined_annotated) = combined_annotated$BEEBASE
combined_annotated = combined_annotated[, !(names(combined_annotated) %in% c("gid", "BEEBASE"))]

write.table(combined_annotated, file = paste("Hb_challenge_drone_ftzf1_gbb", Sys.Date(), sep = "."), row.names = T, col.names = T, quote = F)