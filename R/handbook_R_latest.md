## IO

```R
# IO for table
df <- read.table('path/to/df.txt',header=T,row.names=1,sep='\t')
write.table(df,'path/to/df.txt',sep='\t',col.names=NA)

# IO for sparse matrix

spatial_count <- Matrix::readMM('inputs/CID4535_spatial_count_data/matrix.mtx')
rownames(spatial_count) <- read.table('inputs/CID4535_spatial_count_data/genes.tsv',sep='\t',header=F)$V1
colnames(spatial_count) <- read.table('inputs/CID4535_spatial_count_data/barcodes.tsv',sep='\t',header=F)$V1

writeMM(counts,file='lung_MNP_VERSE_count.mtx')
write.table(rownames(ref_lung@meta.data),file='lung_MNP_VERSE_barcodes.tsv',sep='\t',row.names = F,col.names=F,quote = F)

# RDS
readRDS('file.rds')
saveRDS(obj,'file.rds')

```

## Understand shortcut in Rstudio

* To comment out, using `shift + command + C`

* To clear console, using `ctrl + L`

* To run selection, using `command + enter`.

* clear variable

```R
rm(var1)
rm(list=ls())
```


## Dataframe operation

```R
### Selection can be using (1) index (2) name (3) boolean, negative selection is just add "-" before c
a_data_frame[,c('x','z')]
a_data_frame[,c(1,2)]
a_data_frame[,c(True,False)]  

### access
df$column1
df[['column1']]
```


## baseR functions

```R
table(df$column1) # return a table class
seq(1,10,1)
```



