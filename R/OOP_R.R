# S3 object
s <- list(name="DataFlair",age=29,GPA=4.0)
class(s) <- "student"

s$age  # how to access S3 object's attribute

getGPA <- function(obj){
  UseMethod("getGPA")
}

getGPA.default <- function(obj){
  cat("This is a generic function\n")
}

getGPA.student <- function(obj){
  cat("Total GPA is", obj$GPA,"\n")
}

getGPA(s)

# S4 object

Agent <- setClass(
  "Agent",
  slots = c(location = "numeric", velocity = "numeric", active = "logical"),
  prototype = list(location = c(0.0,0.0),active = T, velocity = c(100.0,0.0)))

a <- new("Agent",location=c(1,5),active = F, velocity = c(9,0))
is.object(a)
isS4(a)
slotNames(a)
slot(a,"location") <- c(1,5)

a@location  # how to access S4 object's attribute

setGeneric("getLocation",function(object){
  standardGeneric("getLocation")
})

setMethod("getLocation",signature(object="Agent"),function(object){
  object@location
})

showMethods("getLocation")
showMethods(class="Agent")

getLocation(a)

# don't care inheritance for now