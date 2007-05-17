setMethod("geneSetType",
          signature=signature(object="GeneIdentifierType"),
          function(object) object@type)

setMethod("show",
          signature=signature(object="GeneIdentifierType"),
          function(object) cat("geneSetType:", geneSetType(object), "\n"))
