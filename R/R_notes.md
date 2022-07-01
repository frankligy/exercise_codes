## Understand where R and its dependencies are installed

```R
R.home()
# /Library/Frameworks/R.framework/Resources
installedPackages <- installed.packages()
# dependencies are either:
## (a) /Library/Frameworks/R.framework/Versions/4.2/Resources/library   [baseR]
## (B) /Users/E0532183/Library/R/x86_64/4.2/library [no renv project]
## (C) /Users/E0532183/cell_typing/renv/library/R-4.2/x86_64-apple-darwin17.0 [renv project]
## (d) /opt/anaconda3/envs/cytospace/lib/R/library [conda] 
sessionInfo()
```

## Understand shortcut in Rstudio

* To comment out, using `shift + command + C`

* To clear console, using `ctrl + L`

* To run selection, using `command + enter`.


## Install R packages

```R
install.packages('devtools')
install.packages("BiocManager")
install.packages('Seurat')
BiocManager::install("SummarizedExperiment")
devtools::install_github('YingMa0107/CARD')
```

