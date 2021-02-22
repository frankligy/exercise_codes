"This script serve as a cheatsheat for base R function"

# make sure there's no dash, punctuation in R col or row name, dot is preferred

# a matrix times a vector can be very tricky, make sure you compute them column-wise, refer
# to TPM functions

# if you don't need name of a vector
a <- c(1,2,3)
names(a) <- c('i1','i2','i3')
a <- unname(a)


# rep command
a <- rep(1,5)

# sample command
sample(c(1:5),3,replace = T)

# cbind function is useful to add a column and combat write table leading empty issue
df <- data.frame(a=c(1,2,3),b=c(5,6,7))
rownames(df) <- c('r1','r2','r3')
df <- cbind("my_column"=c(7,8,9),df)

# write.table, if only need a leading space, set col.names = NA
write.table(df,"/Users/ligk2e/Desktop/test.txt",sep='\t',quote=F,row.names = F)

# rowSums, colSums
rowSums(df,dim=1)
colSums(df,dim=1)

# base plot
boxplot(df)
hist(df$a)
plot(c(1:3),df$b)

# sapply, lapply, apply, tapply

# sapply return a vector


# string in R
# strsplit function

# data frame function in base R
# duplicated


# generate random number
# rnorm, runif


