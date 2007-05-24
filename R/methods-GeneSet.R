.constructors("GeneSet",
              required=c("setName", "setIdentifier"))

setMethod("initialize",
          signature=signature(.Object="GeneSet"),
          function(.Object, .Template=.Object, ...,
                   ## additional args, manipulated by method
                   setIdentifier=.Template@setIdentifier,
                   setName=.Template@setIdentifier,
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

.SETTERS_GeneSet <- .GETTERS_GeneSet

.setters("GeneSet", .SETTERS_GeneSet)

## Logic operations

.checkGeneSetLogicTypes <- function(x, y, functionName) {
    if (!extends(class(setType(x)), class(setType(y))))
        stop(functionName, " requires identical GeneSet setTypes;",
             "\n\tgot: ", class(setType(x)), ", ", class(setType(y)))
}

.geneSetIntersect <- function(x, y) {
   .checkGeneSetLogicTypes(x, y, "'&' or 'intersect'")
    new(class(x), x,
        setIdentifier=setIdentifier(x),
        setName=.glue(setName(x), setName(y), "&"),
        urls=.unique(urls(x), urls(y)),
        genes=intersect(genes(x), genes(y)))
}

.geneSetUnion <- function(x, y) {
    .checkGeneSetLogicTypes(x, y, "'|' or 'union'")
    new(class(x), x,
        setIdentifier=setIdentifier(x),
        setName=.glue(setName(x), setName(y), "|"),
        urls = .unique(urls(x), urls(y)),
        genes=union(genes(x), genes(y)))
}

setMethod("intersect",
          signature=signature(x="GeneSet", y="GeneSet"),
          .geneSetIntersect)

setMethod("union",
          signature=signature(x="GeneSet", y="GeneSet"),
          .geneSetUnion)

setMethod("&",
          signature=signature(e1="GeneSet", e2="GeneSet"),
          function(e1, e2) .geneSetIntersect(e1, e2))

setMethod("&",
          signature=signature(e1="GeneSet", e2="character"),
          function(e1, e2) {
              genes <- intersect(genes(e1), e2)
              new(class(e1), e1,
                  setIdentifier=setIdentifier(e1),
                  setName=.glue(setName(e1), "<character>", "&"),
                  genes=genes)
          })

setMethod("|",
          signature=signature(e1="GeneSet", e2="GeneSet"),
          function(e1, e2) .geneSetUnion(e1, e2))

setMethod("|",
          signature=signature(e1="GeneSet", e2="character"),
          function(e1, e2) {
              genes <- union(genes(e1), e2)
              new(class(e1), e1,
                  setIdentifier=setIdentifier(e1),
                  setName=.glue(setName(e1), "<character>", "|"),
                  genes=genes)
          })

setMethod("Logic",
          signature=signature(e1="character", e2="GeneSet"),
          function(e1, e2) callGeneric(e2, e1))

setMethod("setdiff",
          signature=signature(x="GeneSet", y="GeneSet"),
          function(x, y) {
              .checkSetTypes(x, y, "'setdiff'")
              genes=setdiff(genes(x), genes(y))
              new(class(x), x,
                  setIdentifier=setIdentifier(x),
                  setName=.glue(setName(x), setName(y), "-"),
                  genes=setdiff(genes(x), genes(y)),
                  creationDate=date())
          })

## show

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
