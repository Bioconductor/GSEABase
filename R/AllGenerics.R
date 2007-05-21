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
