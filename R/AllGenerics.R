## many generics are constructed in code; see AAA.R and methods-* 

## GeneIdentifierType

setGeneric("mapIdentifiers",
           function(from, to, what)
           standardGeneric("mapIdentifiers"))

setGeneric("validIdentifiers",
           function(identifier, genes) {
               standardGeneric("validIdentifiers")
           })

## GeneSet

setGeneric("setType<-",
           function(object, value) standardGeneric("setType<-"))

## GeneColorSet

setGeneric("coloring",
           function(object, ...) standardGeneric("coloring"))

setGeneric("coloring<-",
           function(object, ..., value) standardGeneric("coloring<-"))

##

setGeneric("intersect",
           function(x, y) standardGeneric("intersect"))

setGeneric("union",
           function(x, y) standardGeneric("union"))

setGeneric("setdiff",
           function(x, y) standardGeneric("setdiff"))

setGeneric("Logic")
