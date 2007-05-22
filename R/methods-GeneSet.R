.constructors("GeneSet")

setMethod("initialize",
          signature=signature(.Object="GeneSet"),
          function(.Object, .Template=.Object, ...,
                   ## required, even when 'cloning'
                   setIdentifier, setName,
                   ## additional args, manipulated by method
                   shortDescription=.Template@shortDescription,
                   longDescription=.Template@longDescription,
                   organism=.Template@organism,
                   creationDate=date()) {
              callNextMethod(.Object, .Template, ...,
                             setIdentifier=mkScalar(setIdentifier),
                             setName=mkScalar(setName),
                             shortDescription=mkScalar(shortDescription),
                             longDescription=mkScalar(longDescription),
                             organism=mkScalar(organism),
                             creationDate = creationDate)
          })

.GETTERS_GeneSet <- c(setType="type", "genes", "setIdentifier",
                      setName="setName", description="shortDescription",
                      "longDescription", "organism", "pubMedIds", "urls",
                      "contributor", setVersion="version",
                      "creationDate", "collectionType")

.getters("GeneSet", .GETTERS_GeneSet)

## other methods

setMethod("show",
          signature=signature(object="GeneSet"),
          function(object) {
              cat("setName: ", setName(object), "\n",
                  "setIdentifier: ", setIdentifier(object), "\n", sep="")
              cat("genes:",
                  paste(selectSome(genes(object), maxToShow=4), collapse=", "),
                  paste("(total: ", length(genes(object)), ")\n",
                        sep=""),
                  sep=" ")
              show(setType(object))
              show(collectionType(object))
              cat("description: ", description(object), "\n",
                  if(nchar(longDescription(object))!=0 &&
                     longDescription(object) !=  description(object)) {
                      "  (longDescription available)\n"
                  },
                  "organism: ", organism(object), "\n",
                  "pubMedIds: ", pubMedIds(object), "\n",
                  "urls: ", paste(selectSome(urls(object), maxToShow=3),
                                  collapse="\n      "), "\n",
                  "contributor: ", contributor(object), "\n",
                  "setVersion: ",
                  sep="")
              show(setVersion(object))
              cat("creationDate: ", creationDate(object), "\n", sep="")
          })

