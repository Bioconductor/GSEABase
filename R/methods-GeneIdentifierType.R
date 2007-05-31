## GeneIdentifierType

.CONSTRUCTORS_GeneIdentifierType <- 
    names(getSubclasses(getClass("GeneIdentifierType")))

.constructors_Simple(.CONSTRUCTORS_GeneIdentifierType[.CONSTRUCTORS_GeneIdentifierType!="AnnotationIdentifier"])

.constructors_Simple("AnnotationIdentifier", required="annotation")

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

.annotationMapper <- function(from, tag, what) {
    genes <- genes(what)
    pkg <- annotation(from)
    map <- paste(pkg, tag, sep="")
    if (!suppressWarnings(.requireQ(pkg)))
        stop(sprintf("cannot load annotation package '%s'",
                     pkg))
    if (!exists(map))
        stop(sprintf("map '%s' not found in annotation package '%s'",
                     map, pkg))
    ngenes <- mget(genes, get(map))
    if (!all(sapply(ngenes, length)==1))
        stop(sprintf("map '%s' is not 1:1 in annotation package '%s'",
                     map, pkg))
    as.character(unlist(ngenes))
}

## construct these programatically?

setMethod("mapIdentifiers",
          signature=signature(
            from="NullIdentifier", to="GeneIdentifierType",
            what="GeneSet"),
          function(from, to, what) {
              new(class(what), what, type=to)
          })
          

setMethod("mapIdentifiers",
          signature=signature(
            from="AnnotationIdentifier", to="GeneIdentifierType",
            what="GeneSet"),
          function(from, to, what) {
              new(class(what), what,
                  genes=.annotationMapper(from, what,
                    toupper(setType(to))),
                  type=to)
          })

setMethod("mapIdentifiers",
          signature=signature(
            from="AnnotationIdentifier", to="character",
            what="GeneSet"),
          function(from, to, what) {
              type <- tryCatch({
                  res <- do.call(to, list())
                  if (!is(res, "AnnotationIdentifier"))
                      NullIdentifier()
                  else
                      res
              }, error = function(err) NullIdentifier())
              new(class(what), what,
                  genes=.annotationMapper(from, to, what),
                  type=type)
          })

setMethod("mapIdentifiers",
          signature=signature(
            from="AnnotationIdentifier",
            to="EntrezIdentifier",
            what="GeneSet"),
          function(from, to, what) {
              new(class(what), what,
                  genes=.annotationMapper(from, "ENTREZID", what),
                  type=to)
          })
              
setMethod("show",
          signature=signature(object="AnnotationIdentifier"),
          function(object) {
              cat("setType:", setType(object),
                  if (nchar(annotation(object))>0) {
                      paste("(", annotation(object), ")", sep="")
                  }, "\n")
          })
