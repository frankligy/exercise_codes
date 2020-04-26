#!/usr/local/bin/Rscript
library(DoubletDecon,lib.loc = "/Library/Frameworks/R.framework/Versions/3.6/Resources/library")
args<- commandArgs(trailingOnly=TRUE)
filename <- args[1]
removeCC <- as.logical(as.numeric(args[2]))
rhop <- as.numeric(args[3])
centroids <- as.logical(as.numeric(args[4]))
num_doubs <- as.numeric(args[5])
only50 <- as.logical(as.numeric(args[6]))
min_uniq <- as.numeric(args[7])

cat(c(filename,removeCC,rhop,centroid,num_doubs,only50,min_uniq))

Main_Doublet_Decon("/Users/ligk2e/Downloads/DoubletDecon-master/ICGS_expression.txt",
                   "/Users/ligk2e/Downloads/DoubletDecon-master/ICGS_groups.txt",
                   filename,
                   "/Users/ligk2e/Downloads/DoubletDecon-master/output/",
                   fullDataFile=NULL,
                   removeCC=removeCC,
                   species="hsa",
                   rhop=rhop,
                   write=TRUE,
                   PMF=TRUE,
                   useFull=FALSE,
                   heatmap=FALSE,
                   centroids=centroids,
                   num_doubs=num_doubs,
                   only50=only50,
                   min_uniq=min_uniq,
                   nCores=-1)


