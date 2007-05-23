setGeneric("geneSetType",
           function(object) standardGeneric("geneSetType"))

setGeneric("collectionType",
           function(object) standardGeneric("collectionType"))

setGeneric("GeneSet",
           signature=c("type"),
           function(type, ...) standardGeneric("GeneSet"))

setGeneric("GeneColorSet",
           signature=c("type"),
           function(type, ...) standardGeneric("GeneColorSet"))

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
