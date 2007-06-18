## many generics are constructed in code; see AAA.R and methods-* 

## GeneIdentifierType

setGeneric("mapIdentifiers",
           function(what, to, from, ...)
           standardGeneric("mapIdentifiers"),
           useAsDefault=function(what, to, from, ...) {
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

setGeneric("geneIdType<-",
           function(object, value) standardGeneric("geneIdType<-"))

setGeneric("intersect",
           function(x, y) standardGeneric("intersect"))

setGeneric("union",
           function(x, y) standardGeneric("union"))

setGeneric("setdiff",
           function(x, y) standardGeneric("setdiff"))

setGeneric("Logic")

## GeneColorSet

setGeneric("coloring",
           function(object, ...) standardGeneric("coloring"))

setGeneric("coloring<-",
           function(object, ..., value) standardGeneric("coloring<-"))

## GeneSetCollection

setGeneric("GeneSetCollection",
           function(object, ...) standardGeneric("GeneSetCollection"))
