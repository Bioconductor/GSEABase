.constructors("GeneColorSet", required=c("phenotype"))

setMethod("initialize",
          signature=signature(.Object="GeneColorSet"),
          function(.Object, .Template=.Object, ...,
                   ## additional args
                   genes=.Template@genes,
                   phenotype=.Template@phenotype,
                   geneColor=factor(character(length(genes))),
                   phenotypeColor=factor(character(length(genes)))) {
              callNextMethod(.Object, .Template, ...,
                             genes=genes,
                             phenotype=mkScalar(phenotype),
                             geneColor=geneColor,
                             phenotypeColor=phenotypeColor)
          })

.GETTERS_GeneColorSet <- c("phenotype", "geneColor", "phenotypeColor")

.getters("GeneColorSet", .GETTERS_GeneColorSet)

.SETTERS_GeneColorSet <- .GETTERS_GeneColorSet

.setters("GeneColorSet", .SETTERS_GeneColorSet)

setAs("GeneSet", "GeneColorSet",
      function(from) {
          new("GeneColorSet",
              ## clone template
              from,
              ## required args
              setName=setName(from),
              setIdentifier=setIdentifier(from),
              phenotype="undefined")
      })

setMethod("coloring",
          signature=signature(object="GeneColorSet"),
          function(object, ...) {
              data.frame(geneColor=geneColor(object),
                         phenotypeColor=phenotypeColor(object),
                         row.names=genes(object))
          })

setReplaceMethod("coloring",
                 signature=signature(
                   object="GeneColorSet",
                   value="data.frame"),
                 function(object, ..., value) {
                     if (!all(row.names(value) %in% genes(object)))
                         stop("'data.frame' row.names must all be gene symbols")
                     if (nrow(value) != length(genes(object)))
                         stop("'data.frame' must define colors for all genes")
                     if (length(colnames(value)) !=2 ||
                                !all(c("geneColor", "phenotypeColor") %in%
                                     colnames(value)))
                         stop("'data.frame' must only 'geneColor' and 'phenotypeColor' columns")
                     new(class(object), object,
                         genes=row.names(value),
                         geneColor=value[["geneColor"]],
                         phenotypeColor=value[["phenotypeColor"]],
                         setName=setName(object),
                         setIdentifier=setIdentifier(object),
                         phenotype=phenotype(object))
                 })

## other methods

setMethod("show",
          signature=signature(object="GeneColorSet"),
          function(object) {
              callNextMethod()
              cat("phenotype:", phenotype(object), "\n")
              cat("geneColor: ",
                  paste(selectSome(as.character(geneColor(object)),
                                   maxToShow=4),
                        collapse=", "),
                  "\n  levels: ", paste(levels(geneColor(object)),
                                        collapse=", "), "\n",
                  "phenotypeColor: ",
                  paste(selectSome(as.character(phenotypeColor(object)),
                                   maxToShow=4),
                        collapse=", "),
                  "\n  levels: ", paste(levels(phenotypeColor(object)),
                                        collapse=", "), "\n",
                  sep="")
          })
