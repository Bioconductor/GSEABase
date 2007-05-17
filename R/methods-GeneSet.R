.constructors("GeneSet")

.getters("GeneSet",
         c(geneSetType = "type", "genes", "uniqueIdentifier",
         geneSetName = "setName", "shortDescription",
         "longDescription", "organism", "pubMedIds", "urls",
         "contributor", geneSetVersion = "version", "creationDate",
         "collectionType"))

## other methods

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
