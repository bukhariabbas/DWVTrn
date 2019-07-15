# This function does enrichments of TFs targets in gene lists. 
TRN_enrichment = function(files, TRN, universe_size, tfs_cutoff, output_path)
{
  stopifnot(require(reshape))
  stopifnot(require(dplyr))
  #setting the path for reading files and also for outputting files.
  setwd(output_path)
  
  #Converting TRN to TFs -> targets lists
  TRN = aggregate(targets ~ tfs, TRN, toString)
  
  #Triming any character or spaces in target names.
  trim = function (x) gsub("^\\s+|\\s+$", "", x)
  
  # this functions computes enrichment of each TFs in the provide Genelist. Also this function creates a sig_targets 
  # output directory to print overlaping genes between TFs targets and Genelist.
  # After getting an enrichment of TFs in each gene list this then adjusts for multiple testing correcting using bonferroni method.
  getTRN_enrichment = function(TRN, GeneList, universe, filename)
  {
    p = c()
    TF_names = c()
    j = 1
    tf_sig_targets = c()
    for(i in 1:nrow(TRN))
    {
      sample1 = as.character(GeneList)
      sample2 = strsplit(TRN$targets[i], ",")[[1]]
      
      if(length(sample2) >= cutoff) # sample2 contains targets per TF
      {
        sample2 = trim(sample2)
        sample2 = unique(sample2)
        overlap = length(intersect(sample1, sample2))
        tf_sig_targets[j] = toString(intersect(sample1, sample2))
        p[j] = phyper(overlap, length(sample1), universe-length(sample1), length(sample2), lower.tail=F)
        #lower.tail: logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
        TF_names[j] = as.character(TRN$tfs[i])
        j = j + 1
      }
      
    }
    head(tf_sig_targets)
    names(tf_sig_targets) = TF_names
    tf_sig_targets = as.data.frame(tf_sig_targets)
    ### Creat Results Directory.
    if(dir.exists(file.path(paste(getwd(),"/","sig_targets","_",Sys.Date(), sep = ""))) == FALSE)
    {
      dir.create(file.path(paste(getwd(),"/","sig_targets","_",Sys.Date(), sep = "")) ,showWarnings = TRUE)
      setwd(file.path(paste(getwd(),"/","sig_targets","_",Sys.Date(), sep = "")))
    } else {
      setwd(file.path(paste(getwd(),"/","sig_targets","_",Sys.Date(), sep = "")))
    }
    write.table(tf_sig_targets, file = filename, sep = "\t", col.names = F, row.names = T, quote = F)
    p = p.adjust(p, method = "bonferroni")
    p = as.data.frame(p)
    
    p$tfs = TF_names
    colnames(p) = c(filename, "tfs")
    
    return(p)
  }
  
  # This Pvalues object holds significant enrichments.
  Pvalues = c()
  # This for loop is the interface which uses enrichment function iteratively for each gene list provided.
  for( i in 1: length(files))
  {
    setwd(path)
    filename = strsplit(files[i], "/")[[1]]
    filename = filename[length(filename)]
    iFile = read.delim(files[i], header = F)
    #  iFile = iFile[which(iFile$CFDR < 0.1), c("ensembl_gene_id")]
    if(length(iFile$V1) >= cutoff)
    {
      if (i == 1)
      {
        Pvalues = getTRN_enrichment(TRN, iFile$V1, universe, filename)
      }
      else
      {
        Pvalues = plyr::join(Pvalues, getTRN_enrichment(TRN, iFile$V1, universe, filename), by = "tfs")
      }
    }
  }
  
  rownames(Pvalues) = Pvalues$tfs
  Pvalues = Pvalues[,c(-2)]
  return(Pvalues)
}
