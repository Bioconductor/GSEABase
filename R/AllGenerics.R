setGeneric("geneSetType",
           function(object) standardGeneric("geneSetType"))

setGeneric("collectionType",
           function(object) standardGeneric("collectionType"))

setGeneric("GeneSet",
           function(type, ..., creationDate = date()) {
               standardGeneric("GeneSet")
           }, signature=c("type"))
