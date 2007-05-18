.constructors("GeneSet")

setMethod("initialize", "GeneSet",
          function(.Object, ...,
                   setIdentifier, setName,
                   shortDescription="", longDescription=shortDescription,
                   organism="",
                   creationDate=date()) {
              callNextMethod(.Object, ...,
                             setIdentifier=mkScalar(setIdentifier),
                             setName=mkScalar(setName),
                             shortDescription=mkScalar(shortDescription),
                             longDescription=mkScalar(longDescription),
                             organism=mkScalar(organism),
                             creationDate = creationDate)
          })

.GETTERS_GeneSet <- c(geneSetType="type", "genes", "setIdentifier",
                      geneSetName="setName", "shortDescription",
                      "longDescription", "organism", "pubMedIds", "urls",
                      "contributor", geneSetVersion="version",
                      "creationDate", "collectionType")

.getters("GeneSet", .GETTERS_GeneSet)

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
                  "setIdentifier: ", setIdentifier(object), "\n",
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
