.constructors("GeneColorSet")

setMethod("initialize",
          signature=signature(.Object="GeneColorSet"),
          function(.Object, .Template=.Object, ...,
                   ## required arg, even when cloning
                   phenotype,
                   ## additional args
                   genes=.Template@genes,
                   geneColor=factor(character(length(genes))),
                   phenotypeColor=factor(character(length(genes)))) {
              callNextMethod(.Object, .Template, ...,
                             genes=genes,
                             phenotype=mkScalar(phenotype),
                             geneColor=geneColor,
                             phenotypeColor=phenotypeColor)
          })

.getters("GeneColorSet",
         c("phenotype", "geneColor", "phenotypeColor"))

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

## other methods

setMethod("show",
          signature=signature(object="GeneColorSet"),
          function(object) {
              callNextMethod()
              cat("phenotype:", phenotype(object), "\n")
              cat("geneColor:",
                  paste(selectSome(geneColor(object)), collapse=", "),
                  "\n  levels:", levels(geneColor(object)), "\n",
                  "phenotypeColor:",
                  paste(selectSome(phenotypeColor(object)), collapse=", "),
                  "\n  levels:", levels(phenotypeColor(object)), "\n")
          })
