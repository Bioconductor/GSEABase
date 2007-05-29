## GeneIdentifierType

.CONSTRUCTORS_GeneIdentifierType <- 
    names(getSubclasses(getClass("GeneIdentifierType")))

.constructors_Simple(.CONSTRUCTORS_GeneIdentifierType)

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

## AnnotationIdentifier

setMethod("initialize",
          signature=signature(.Object="AnnotationIdentifier"),
          function(.Object, .Template=.Object, ...,
                   annotation=.Template@annotation) {
              callNextMethod(.Object, .Template, ...,
                             annotation = mkScalar(annotation))
          })

.GETTERS_AnnotationIdentifier <- c("annotation")

.getters("AnnotationIdentifier", .GETTERS_AnnotationIdentifier)

setMethod("mapIdentifiers",
          signature=signature(
            from="AnnotationIdentifier",
            to="EntrezIdentifier",
            what="GeneSet"),
          function(from, to, what) {
              genes <- genes(what)
              pkg <- annotation(from)
              genemap <- paste(pkg, "ENTREZID", sep="")
              if (!suppressWarnings(.requireQ(pkg)))
                  stop(sprintf("cannot load annotation package '%s'",
                               pkg))
              if (!exists(paste(pkg, "ENTREZID", sep="")))
                  stop(sprintf("cannot find map '%s' in annotation package '%s'",
                               genemap, pkg))
              ngenes <- mget(genes, get(genemap))
              if (!all(sapply(ngenes, length)==1))
                  stop(sprintf("mapping between '%s' and '%s' is not 1:1 in '%s'",
                               "AnnotationIdentifier",
                               "EntrezIdentifier", genemap))
              new(class(what), what,
                  genes=as.character(unlist(ngenes)), type=to)
          })

setMethod("show",
          signature=signature(object="AnnotationIdentifier"),
          function(object) {
              cat("setType:", setType(object),
                  if (nchar(annotation(object))>0) {
                      paste("(", annotation(object), ")", sep="")
                  }, "\n")
          })
