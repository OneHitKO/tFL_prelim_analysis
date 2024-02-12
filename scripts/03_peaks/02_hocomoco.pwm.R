library(universalmotif)
library(TFBSTools)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# always set seed!!!
set.seed(105)


# Get PWM of TF motifs for downstream peak anno; use HOCOMOCO v12
# import files and jsonl
dir = "/g/zaugg/zaugg_shared/databases/TFs/motifs/HOCOMOCO/v12/"
pwmfiles = list.files(path = paste0(dir,"H12CORE_pwm/pwm/"), full.names = T)
pwmnames = gsub("\\.pwm","",basename(pwmfiles))

# separate the score from TF name
pwmnames = gsub("\\.H12CORE\\.","_",pwmnames)

# read annotation from json file
lines = readLines(paste0(dir,"H12CORE_annotation.jsonl"))
lines = lapply(json, jsonlite::fromJSON)
lines = lapply(lines, unlist)
anno = bind_rows(lines)

# import pwm
pwmlist = map(
  pwmfiles,
  ~ read_matrix(
    .x, 
    type = "PWM",
    positions = "rows"
  )
)

names(pwmlist) = pwmnames

# gsub name
pwmlist = lapply(pwmlist, function(x){
  new = gsub(">","",x@name)
  x@name = new
  return(x)
})

# convert 
finalPWM = convert_motifs(pwmlist, class = "TFBSTools-PWMatrix")

# add ID
finalPWM = lapply(finalPWM, function(x){
  x@ID = x@name
  return(x)
})

finalPWM = do.call(PWMatrixList,finalPWM)

# save
qs::qsave(finalPWM,"rdata/03_02_finalPWM.qs")

