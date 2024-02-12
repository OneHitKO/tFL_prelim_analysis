#!/bin/bash
#
# ================================================= #
# numbat: prepare allele data
# ================================================= #
#
#SBATCH --job-name=numbat
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kristy.ou@embl.de
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=99:00:00
#SBATCH --error=.slurmR/%j_phase.err
#SBATCH --output=.slurmR/%j_phase.out

# This script counts alleles and phase SNPs for each sample.
# Run from project dir with sbatch scripts/04_mutation/01_prepare.numbat.sh

module purge
module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2

# activate conda environment
source /home/kou/.bashrc
mamba activate numbat

# define input dir
indir="/g/saka/Kristy/projects/composite"

# define LN
LN="LN0025 LN0027 LN0177 LN0193 LN0438"

# define paths
scriptPath="/g/saka/Kristy/R-lib/4.2.2-foss-2022b/numbat/bin/pileup_and_phase.R"
bamDir="${indir}/cellranger"
whitelistDir="${indir}/analysis/other/whitelist"
eagleBin="/g/saka/Kristy/software/envs/numbat/bin/eagle"
gmap="/scratch/kou/ref/eagle_genetic_map_hg38_withX.txt.gz"
snpVCF="/scratch/kou/ref/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
phasingPanel="/scratch/kou/ref/1000G_hg38"
outdir="/scratch/kou/composite/numbat"

# perform function
for ln in $LN
do
Rscript $scriptPath \
 --label $ln \
 --samples $ln \
 --bams ${bamDir}/${ln}/outs/gex_possorted_bam.bam \
 --barcodes ${whitelistDir}/${ln}_bamCB.tsv \
 --gmap $gmap \
 --eagle $eagleBin \
 --snpvcf $snpVCF \
 --paneldir $phasingPanel \
 --outdir $outdir/${ln} \
 --ncores 64 &
done
wait
 
