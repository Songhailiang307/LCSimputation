#!~/.conda/envs/myrenv/bin/ Rscript
library("STITCH")

Args <- commandArgs(trailingOnly=TRUE)

CHR=as.numeric(Args[1])
BAMLIST=Args[2]
POS=Args[3]




STITCH(chr= CHR, bamlist=BAMLIST, posfile=POS, outputdir="./", K=10, nGen=100, nCores=20,tempdir="/data/songhl/temp")
