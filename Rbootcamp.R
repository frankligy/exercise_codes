#### configurations
# following four lines are borrowed from https://github.com/EDePasquale/DoubletDecon/blob/master/R/Main_Doublet_Decon.R
# just for the use of learning
# nice packages: stringr, tidyr, dplyr,foreach,doParallel,plyr,tidyverse(a set of packages), installr, data.table, R.utils, gplots, ggplot2
#              reshape2
options(install.packages.compile.from.source = "never")
if(!("BiocManager" %in% installed.packages()[, "Package"])) install.packages("BiocManager")
library(BiocManager)
if(!("DeconRNASeq" %in% installed.packages()[, "Package"])) BiocManager::install("DeconRNASeq", update=F)
if(!("gplots" %in% installed.packages()[, "Package"])) install.packages("gplots")

# alternative way: using require()
# devtools to intall source code from github
if(!require(devtools)){
    install.packages("devtools")
}
devtools::install_github('EDePasquale/DoubletDecon')   # DoubletDecon is the name of the package

# alternative way: use biocLite (Seems from R > 3.5, we should deprecate this, instead, use BiocManager shown as above)
source("https://bioconductor.org/biocLite.R")
biocLite(c("DeconRNASeq", "clusterProfiler", "hopach", "mygene", "tidyr", "R.utils", "foreach", "doParallel", "stringr"))
install.packages("MCL")

# if you encounter say Warning: package is not available for certain R version, specify your mirror to following
# if still not work, well, might your current version indeed not compatible with the package.
# or, we should use biocManager or devtool to install, just check each packages requirement
install.packages('shinyDirectoryInput',repos='http://cran.us.r-project.org')

# all the installed packages
a = installed.packages() # a will be a matrix with all installed packages and their attributes
b = library()   # b will be a S3 object
b1 = library()$result  # this attribute will ba a matrix, recording the all the available path in searched path
path = .libPaths()  # this is the path that R will search for

# check the version of R and Rstudio
"Rstudio: open RStudio in the upper left corner and see \"About Rstudio\"
 R: in terminal type: R --version  "

# learning R:
"1. base functionality
 2. statistics, it is about math, not about R
 3. visualization, please refer to http://www.sthda.com/english/wiki/data-visualization"


## set of inspect function
class()
typeof()
is.object()
isS4()

# if matrix or data.frame
dim()
colnames()
rownames()

# if list
length()

# R script, take arguments
# Rscript whatever.R 1 2 3 4 5 6 7
# args[1] will be "1"
args<- commandArgs(trailingOnly=TRUE)
filename <- args[1]
removeCC <- as.logical(as.numeric(args[2]))
rhop <- as.numeric(args[3])
centroid <- as.logical(as.numeric(args[4]))
num_doubs <- as.numeric(args[5])
only50 <- as.logical(as.numeric(args[6]))
min_uniq <- as.numeric(args[7])

# RDS file is to sava a single object
save(pbmc,file="")
pbmc = load("")

# data.frame operation
# example1:
# dplyr packages: group_by,select,filter,arrange,mutate,pipe symbol
genes=seuratObject.markers %>% group_by(cluster) %>% top_n(n = num_genes, wt = avg_logFC)

# example2: reorder row-wise and column-wise
genes2=genes[order(genes$cluster),]   # row-wise
expression=expression[row.names(clusters2)]  # column-wise
expression=expression[match(geneOrder, row.names(expression)),]  # match function will return the index that each item in vector1 occur in vector2

# factor is a data structure,using class function to view, labels are the name, it could be converted to data.frame, labels will be rownames
directions <- c("North","East","South","West")
a = factor(directions, levels= c("North", "East", "South", "West"), labels=c("N", "E", "S", "W"))
a1 = as.data.frame(a)
b = factor(c(1,2,3,4))
b1 = as.character(b)
b2 = as.numeric(b1)


























