## many generics are constructed in code; see AAA.R and methods-* 

## GeneIdentifierType

setGeneric("mapIdentifiers",
           function(what, to, from, ..., verbose=FALSE)
           standardGeneric("mapIdentifiers"),
           signature=c("what", "to", "from"),
           useAsDefault=function(what, to, from, ..., verbose=FALSE) {
               if (geneIdType(from) == geneIdType(to)) {
                   warning(sprintf("map from '%s' to '%s': identical types",
                                   geneIdType(from), geneIdType(to)))
                   what
               } else {
                   stop(sprintf("cannot map from '%s' to '%s' on object of class '%s'",
                                geneIdType(from), geneIdType(to),
                                class(what)))
               }
           })

## GeneSet

setGeneric("details",
           function(object) standardGeneric("details"))

setGeneric("geneIdType<-",
           function(object, verbose=FALSE, value)
           standardGeneric("geneIdType<-"),
           signature=c("object", "value"))

setGeneric("intersect",
           function(x, y) standardGeneric("intersect"))

setGeneric("union",
           function(x, y) standardGeneric("union"))

setGeneric("setdiff",
           function(x, y) standardGeneric("setdiff"))

setGeneric("incidence",
           function(x, ...) standardGeneric("incidence"))

## GeneColorSet

setGeneric("coloring",
           function(object, ...) standardGeneric("coloring"))

setGeneric("coloring<-",
           function(object, ..., value) standardGeneric("coloring<-"))

## GeneSetCollection

setGeneric("GeneSetCollection",
           function(object, ..., idType, setType)
           standardGeneric("GeneSetCollection"))
