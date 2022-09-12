# R doesn't have scalar type, for example, in python, we have int type(object), float type(object), which
# is scalar type. In R, everything is a vector, vector is the smallest unit. 

# R datatype
"
vector: can be created using c() function or colon, it can have names
  numeric vector:
    integer vector
    double vector
    complex vector
  character vector
  raw vector: rarely used
  
To index a vector, passing a vector, either numeric, or complementary numeric, or boolean, or names

factor: factor is a vector with label
list: a vector, but each element can be different type of different length
matrix
array
dataframe: each column will be a factor

to check:
class(variable)

since list can contain sublist, list is recursive type, others are atomic type

useful functions:
object.size()
rm(list=ls())

Enviroment is to store variables, it has function-specific environment, and the global enviroment 
Envionment itself is a list-like variable as well.
"

v1 <- c(1,2,3)
names(v1) <- c("first","second","third")
v1[c(1,2)]
v1[c(-3)]
v1[c(T,T,F)]
v1[c("first","second")]

class(v1)

gender <- factor(c("male","female","female","female","male"))
names(gender) <- c("first","second","third","fourth","fifth")
levels(gender)
nlevels(gender)
as.integer(gender)
str(gender)  # structure of an object
attributes(gender)  # all the attributes of an object, a factor object

# both matrix and array are filled in column-wise, different from python, you can change it by
# setting byrow=TRUE

a_matrix <- matrix(
  1:12,
  nrow=4,   # or dim=c(4,3)
  dimnames=list(
    c("one","two","three","four"),
    c("ein","zwei","drei")
  )
)

dim(a_matrix)
nrow(a_matrix)
ncol(a_matrix)
length(a_matrix)  # the product of each dimension
rownames(a_matrix)
colnames(a_matrix)
dimnames(a_matrix)  # a list

# you can cbind and rbind matrix as well
# you can transpose matrix by t() function
# you can do inner multiplication %*% and outer multiplication by %o%
# you can inverse matrix using solve() function


three_d_array <- array(
  1:24,
  dim=c(4,3,2),
  dimnames=list(
    c("one","two","three","four"),
    c("ein","zwei","drei"),
    c("un","deux")
  )
)  

dim(three_d_array)
dimnames(three_d_array)
length(three_d_array)



a_list <- list(
  c(1,2,3,4),
  matrix(1:12,nrow=3)
)

names(a_list) = c("item1","item2")
a_list$item1  
a_list[["item1"]]
a_list["item1"]
# you can use vector based methods but will return a new list, using dollor sign or double bracket can
# directly get the content stored in each slot
unlist(a_list)

a_data_frame <- data.frame(
  x=letters[1:5],
  y=rnorm(5),
  z=runif(5)>0.5
)

rownames(a_data_frame)
colnames(a_data_frame)
dimnames(a_data_frame)

# dataframe related useful function

### Selection can be using (1) index (2) name (3) boolean, negative selection is just add "-" before c
a_data_frame[,c('x','z')]
a_data_frame[,c(1,2)]
a_data_frame[,c(True,False)]  


# very useful function in R

"
-as and -is series
rep
seq
paste
which
identical(v1,v2)   # tolerate an infinitesimal number epsilon
NROW,NCOL
all.equal(v1,v2)
any and all
"

rep(1:5,3)
rep(1:5,each=3)
rep(1:5,times=1:5)
rep(1:5,length.out=7)

seq.int(0.1,0.01,-0.01)

v2 <- c(1,2,3)
which( v2 > 2)  # return the location where a logical vector is true

# coerce v2 to a column vector of dimension (3,1), normaly, nrow and ncol won't work on vector
NROW(v2)  
NCOL(v2)


# compare five special variable, enclose each of them using doubel quotation.
"
inf
-inf
NAN    length would be 1, signify mathematical operation can not be made
NA     length would be 1, signify missing value 
NULL   length would be 0, often occur in list, can be used to remove stuff
"







