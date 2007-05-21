.getters("GeneIdentifierType", c(setType="type"))

setMethod("validIdentifiers",
          signature=signature(
            identifier="GeneIdentifierType"),
          function(identifier, genes) {
              TRUE
          })

setMethod("show",
          signature=signature(object="GeneIdentifierType"),
          function(object) cat("setType:", setType(object), "\n"))
