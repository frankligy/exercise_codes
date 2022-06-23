## This is a comprehensive notes for R

# installation
### first check if the pacakge has been installed, and familiar with where the R package goes
installedPackages <- installed.packages()   # /Users/E0532183/Library/R/x86_64/4.2/library or /Library/Frameworks/R.framework/Versions/4.2/Resources/library or /opt/anaconda3/envs/cytospace/lib/R/library if using conda
libraryPaths <- library()$result
searchPath <- .libPaths()
### auxiliary downloading tool
install.packages('devtools')
install.packages("BiocManager")
### download
BiocManager::install("SummarizedExperiment")
devtools::install_github('YingMa0107/CARD')

### debug
# 1. in where R is installed, there's a Makeconf file can be editted (/Library/Frameworks/R.framework/Resources/etc/Makeconf)
#    to configure the path to the gcc, clang, gfortran compiler.
# 2. we can potentially create a file called ~/.R/Makevars to overwrite the setting in Makeconf, no guarantee will work
# 3. install gfortran or other compilers (https://mac.r-project.org/tools/), add it to the path


# we recommend using renv for environment management
### In that way, R packages will be in: /Users/E0532183/cell_typing/renv/library/R-4.2/x86_64-apple-darwin17.0
### If you want to run it outside of Rstudio: follow the instructions here: https://stackoverflow.com/questions/66311017/rscript-to-use-renv-environment


# Basic R data type
### Please refer to https://towardsdatascience.com/learning-r-data-types-e698d23f8179


# R IO
### 1. read.table (https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/read.table
### 2. Matrix package Matrix::readMM('path_to_mtx.mtx'), and then use rownames and colnames function to assign row,col names
### 3. write.table, make sure to set col.names =NA to avoid leading entry empty issue
### 4. RDS, using readRDS or saveRDS

# R commands
rm(a)
rm(list=ls())
R.home()   # check where is R installed
sessionInfo()  # print all dependencies

# keyboard shortcut on Mac
### comment out: shift + command + C
### clear console: ctrl + L
### Run selection: command + enter



# base R function
seq(1,10,1)

