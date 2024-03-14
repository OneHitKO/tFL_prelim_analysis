library(ArchR)
library(Seurat)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tidyverse)
library(SummarizedExperiment)


# Add Peaks to entire wnn object
# import data ----
archr = loadArchRProject("rdata/archr/prelim_anno_peaks/")
wnn = qs::qread("rdata/02_04_integrated.qs", nthreads = 32)
LN_ID = c("LN0025","LN0027","LN0177","LN0193","LN0438")

# create FragmentObject ----
# get frag files
dir = paste0("/g/saka/Kristy/projects/composite/cellranger/",LN_ID,"/outs")
fragmentFiles = list.files(path=dir, full.names = T, pattern = "fragments.tsv.gz$")

# get whitelisted cells per donor
Idents(wnn) = "Sample"

cellsKeep = map(
  LN_ID,
  ~ WhichCells(wnn, idents = .x) %>%
    gsub("LN.*_","",.)
)

# create FragmentObject
fragObj = map2(
  fragmentFiles,
  cellsKeep,
  ~ CreateFragmentObject(
    path = .x,
    cells = .y,
    validate.fragments = F,
    verbose = T
  )
)

# create MotifObject ----
# get peakset GRanges
peakGR = getPeakSet(archr)

# create unique IDs for peaks
peakID = paste0(seqnames(peakGR),":",start(peakGR),"-",end(peakGR))

# get paths to pwm and positions in peak
peakAnno_path = getPeakAnnotation(archr, name = "HOCOv12")
motifPWM = peakAnno_path$motifs
motifPosition = readRDS(peakAnno_path$Positions)
motifMatch = readRDS(peakAnno_path$Matches)

# get matrix with region x motif matches
motifMatchMat = assay(motifMatch)
rownames(motifMatchMat) = peakID #rows need to match peak mat later!!

# create motif obj, motifs are rows!
motifObj = CreateMotifObject(
  data = motifMatchMat,
  pwm = motifPWM,
  motif.names = names(motifPWM),
  positions = motifPosition
)

# create ChromatinAssay with peaks ----
# get peak matrix with peaks vs cell
peakSE = getMatrixFromProject(archr,"PeakMatrix")
peakMat = assay(peakSE)
rownames(peakMat) = peakID
colnames(peakMat) = gsub("#","_",colnames(peakMat))

# get annotation
gencode = rtracklayer::import("/g/zaugg/zaugg_shared/annotations/hg38/Gencode_v45/gencode.v45.annotation.gtf")

# create assay
chromAssay = CreateChromatinAssay(
  counts = peakMat,
  ranges = peakGR,
  fragments = fragObj,
  genome = "hg38",
  annotation = gencode,
  validate.fragments = F,
  verbose = T
)

# add motifs
chromAssay@motifs = motifObj

# Add to wnn seurat object ----
wnn[["Peaks"]] = chromAssay

# add metadata
peakMeta = colData(peakSE)
rownames(peakMeta) = gsub("#","_", rownames(peakMeta))
ord = match(Cells(wnn), rownames(peakMeta))
peakMeta = peakMeta[ord,]

wnn$FRIP = peakMeta$FRIP

# diet then save
diet = DietSeurat(
  wnn,
  assays = c(
    "RNA",
    "SCT",
    "Peaks",
    "GeneScoreImputed"
  ),
  dimreducs = c(
    "rpca.rnaPCA.wnnUMAP",
    "rnaPCA",
    "rpca.rnaPCA",
    "rpca.rnaUMAP",
    "atacUMAP",
    "integratedLSI"
  ),
  graphs = c(
    "rpca.rnaPCA.wknn",
    "rpca.rnaPCA.wsnn"
  )
)

# reset idents
Idents(diet) = "prelim_anno"

qs::qsave(diet, "rdata/03_04_wnn.wPeaks.qs", nthreads = 32)

# Subset B and T ----
# get levels
Bcells = levels(diet)[grep("B_",levels(diet))]
Tcells = levels(diet)[grep("T_",levels(diet))]

# NEED TO REMOVE MOTIF OBJECT to subset! urghhh
# TO REPEAT AFTER DL SIGNAC FROM DEV BRANCH!
diet[["Peaks"]]@motifs = NULL

# subset
diet_B = subset(diet, idents = Bcells)
diet_T = subset(diet, idents = Tcells)

# add chromVAR ----
archr_B = loadArchRProject("rdata/archr/prelim_anno_B")
archr_T = loadArchRProject("rdata/archr/prelim_anno_T")

# get motif deviations and z scores for B
motifSE_B = getMatrixFromProject(archr_B, "HOCOv12Matrix", verbose = T)
devMat_B = assay(motifSE_B, "deviations")
zMat_B = assay(motifSE_B, "z")

# rename columns
colnames(devMat_B) = gsub("#","_",colnames(devMat_B))
colnames(zMat_B) = gsub("#","_",colnames(zMat_B))

# get motif deviations T
motifSE_T = getMatrixFromProject(archr_T, "HOCOv12Matrix", verbose = F)
devMat_T = assay(motifSE_T, "deviations")
zMat_T = assay(motifSE_T, "z")

# rename columns
colnames(devMat_T) = gsub("#","_",colnames(devMat_T))
colnames(zMat_T) = gsub("#","_",colnames(zMat_T))

# create assay
chrVar_B = CreateAssay5Object(
  counts = devMat_B,
  data = zMat_B
)

chrVar_T = CreateAssay5Object(
  counts = devMat_T,
  data = zMat_T
)

# add to object
diet_B[["chromVAR"]] = chrVar_B
diet_T[["chromVAR"]] = chrVar_T

# save for future
qs::qsave(diet_B, "rdata/03_04_wnnB.qs", nthreads = 32)
qs::qsave(diet_T, "rdata/03_04_wnnT.qs", nthreads = 32)
