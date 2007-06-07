## many generics are constructed in code; see AAA.R and methods-* 

## GeneIdentifierType

setGeneric("mapIdentifiers",
           function(what, to, from, ...)
           standardGeneric("mapIdentifiers"),
           useAsDefault=function(what, to, from, ...) {
               if (setType(from) == setType(to)) {
                   warning(sprintf("map from '%s' to '%s': identical types",
                                   setType(from), setType(to)))
                   what
               } else {
                   stop(sprintf("cannot map from '%s' to '%s' on object of class '%s'",
                                setType(from), setType(to),
                                class(what)))
               }
           })

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
