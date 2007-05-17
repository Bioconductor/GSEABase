.getters("GeneIdentifierType", c(geneSetType="type"))

setMethod("validIdentifiers",
          signature=signature(
            identifier="EntrezIdentifier"),
          function(identifier, genes) {
              
          })

setMethod("show",
          signature=signature(object="GeneIdentifierType"),
          function(object) cat("geneSetType:", geneSetType(object), "\n"))
