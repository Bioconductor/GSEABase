## constructors
setMethod("GeneSet",
          signature=signature(type="missing"),
          function(type, ..., creationDate) {
              new("GeneSet", type = new("Untyped"),
                  creationDate = creationDate)
          })

setMethod("GeneSet",
          signature=signature(type="character"),
          function(type, ..., creationDate) {
              new("GeneSet", type = new(type), ...,
                  creationDate=creationDate)
          })

setMethod("GeneSet",
          signature=signature(type="GeneIdentifierType"),
          function(type, ..., creationDate) {
              new("GeneSet", type = type, ...,
                  creationDate = creationDate)
          })

.getters <- function() {
    ## getter name = slot; missing name uses slot for getter name
    slots <-
        c(geneSetType = "type", "genes", "uniqueIdentifier",
          geneSetName = "setName", "shortDescription",
          "longDescription", "organism", "pubMedIds", "urls",
          "contributor", geneSetVersion = "version", "creationDate",
          "collectionType")
    autoName <- nchar(names(slots))==0
    names(slots)[autoName] <- slots[autoName]
    ## standard getters. 'where' default is topenv(parent.frame())
    ## which on package load is the package name space
    for (i in seq(along=slots)) {
        eval(substitute({
            setGeneric(GENERIC,
                       function(object) standardGeneric(GENERIC))
            setMethod(GENERIC,
                      signature=signature(object="GeneSet"),
                      function(object) slot(object, SLOT))
        }, list(GENERIC = names(slots)[[i]],
                SLOT = slots[[i]])))
    }
}

.getters()

setMethod("show",
          signature=signature(object="GeneSet"),
          function(object) {
              cat("genes:", head(genes(object)),
                  paste("(length: ", length(genes(object)), ")\n",
                        sep=""),
                  sep="")
              show(geneSetType(object))
              show(collectionType(object))
              cat(
                  "geneSetName: ", geneSetName(object), "\n",
                  "uniqueIdentifier: ", uniqueIdentifier(object), "\n",
                  "shortDescription: ", shortDescription(object), "\n",
                  "longDescription: ",
                  if(nchar(longDescription(object))==0) "not ",
                  "available", "\n",
                  "organism: ", organism(object), "\n",
                  "pubMedIds: ", pubMedIds(object), "\n",
                  "urls: ", urls(object), "\n",
                  "contributor: ", contributor(object), "\n",
                  "geneSetVersion: ",
                  sep="")
              show(geneSetVersion(object))
              cat("creationDate: ", creationDate(object), "\n", sep="")
          })
