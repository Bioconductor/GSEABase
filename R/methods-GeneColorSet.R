.constructors("GeneColorSet")

setMethod("initialize", "GeneColorSet",
          function(.Object, ...,
                   genes=character(0), phenotype,
                   geneColor=factor(character(length(genes))),
                   phenotypeColor=factor(character(length(genes)))) {
              callNextMethod(.Object, ...,
                             genes=genes,
                             phenotype=mkScalar(phenotype),
                             geneColor=geneColor,
                             phenotypeColor=phenotypeColor)
          })

.getters("GeneColorSet",
         c("phenotype", "geneColor", "phenotypeColor"))

## other methods

setMethod("show",
          signature=signature(object="GeneColorSet"),
          function(object) {
              cat("genes:", head(genes(object)),
                  paste(" (length: ", length(genes(object)), ")\n",
                        sep=""),
                  sep="")
              cat("phenotype:", phenotype(object), "\n")
              cat("geneColor:",
                  head(geneColor(object)),
                  "\n  levels:", levels(geneColor(object)), "\n")
              cat("phenotypeColor:",
                  head(phenotypeColor(object)),
                  "\n  levels:", levels(phenotypeColor(object)), "\n")
              show(geneSetType(object))
              show(collectionType(object))
              cat("geneSetName: ", geneSetName(object), "\n",
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
