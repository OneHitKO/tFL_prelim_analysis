library(tidyverse)
library(Seurat)

# this script creates a list of whitelisted barcodes with or without cell anno depending on software
# notes: CB of bam file includes -1! sometimes this is necessary for software
outdir = "other/whitelist/"

wnn = qs::qread("rdata/02_04_integrated.qs", nthreads = 32)
meta = as_tibble(wnn@meta.data, rownames = "barcodes")
meta$bam_bc = gsub("LN.*_","",meta$barcodes)

# group by sample
meta = meta %>% 
  split(f = as.factor(.$Sample))

# CREATE WHITELIST FOR NUMBAT: barcodes WITH -1 to match bam! 
lapply(seq_along(meta), function(x){
  
  write_tsv(
    as.data.frame(meta[[x]]$bam_bc), 
    paste0(outdir,names(meta)[x],"_bamCB.tsv"), col_names = F, append = F
    )
})

