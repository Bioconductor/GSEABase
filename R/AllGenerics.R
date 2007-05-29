setGeneric("geneSetType",
           function(object) standardGeneric("geneSetType"))

setGeneric("collectionType",
           function(object) standardGeneric("collectionType"))

setGeneric("coloring",
           function(object, ...) standardGeneric("coloring"))

setGeneric("coloring<-",
           function(object, ..., value) standardGeneric("coloring<-"))

setGeneric("validIdentifiers",
           function(identifier, genes) {
               standardGeneric("validIdentifiers")
           })

##

setGeneric("intersect",
           function(x, y) standardGeneric("intersect"))

setGeneric("union",
           function(x, y) standardGeneric("union"))

setGeneric("setdiff",
           function(x, y) standardGeneric("setdiff"))

setGeneric("Logic")
