## GeneIdentifierType

.CONSTRUCTORS_GeneIdentifierType <- 
    names(getSubclasses(getClass("GeneIdentifierType")))

.CONSTRUCTORS_GeneIdentifierType <-
    .CONSTRUCTORS_GeneIdentifierType[!(.CONSTRUCTORS_GeneIdentifierType %in%
                                       c("AnnotationIdentifier"))]

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

.SETTERS_AnnotationIdentifier <-
    .GETTERS_AnnotationIdentifier <- c("annotation")

.getters("AnnotationIdentifier", .GETTERS_AnnotationIdentifier)

.setters("AnnotationIdentifier", .SETTERS_AnnotationIdentifier)

setMethod("show",
          signature=signature(object="AnnotationIdentifier"),
          function(object) {
              cat("setType:", setType(object),
                  if (nchar(annotation(object))>0) {
                      paste("(", annotation(object), ")", sep="")
                  }, "\n")
          })

## mapIdentifiers

setMethod("mapIdentifiers",
          signature=signature(
            what="GeneSet",
            to="GeneIdentifierType",
            from="missing"),
          function(what, to, from, ...) {
              callGeneric(what, to, from=setType(what), ...)
          })

## Null --> X

setMethod("mapIdentifiers",
          signature=signature(
            what="GeneSet",
            to="GeneIdentifierType",
            from="NullIdentifier"),
          function(what, to, from, ...) {
              new(class(what), what, type=to)
          })
          

## AnnotationIdentifier --> X

.checkPackageAndMap <- function(pkg, map) {
    if (!.requireQ(pkg))
        .stopf("cannot load annotation package '%s'", pkg)
    if (!exists(map, where=getNamespace(pkg)))
        .stopf("map '%s' not found in annotation package '%s'",
               map, pkg)
}

.getMappedGenes <- function(genes, mapEnv, map, pkg) {
    ngenes <- mget(genes, mapEnv)
    if (any(length(ngenes) != 1) || any(is.na(ngenes)))
        .warningf("annotation map '%s' is not 1:1 in '%s'",
                  map, pkg)
    ugenes <- unique(unlist(ngenes))
    as.character(ugenes[!is.na(ugenes)])
}

## tag: e.g., ENTREZID
.fromAnnotation <- function(from, tag, what) {
    pkg <- annotation(from)
    map <- paste(gsub("db$", "", pkg), tag, sep="")
    .checkPackageAndMap(pkg, map)
    mapEnv <- get(map, envir=getNamespace(pkg))
    .getMappedGenes(genes(what), mapEnv, map, pkg)
}

setMethod("mapIdentifiers",
          signature=signature(
            what="GeneSet",
            to="GeneIdentifierType",
            from="AnnotationIdentifier"),
          function(what, to, from, ...) {
              tag <- toupper(setType(to))
              new(class(what), what,
                  genes=.fromAnnotation(from, tag, what),
                  type=to)
          })

## X --> AnnotationIdentifier

## tag: e.g., SYMBOL
.toAnnotation <- function(tag, to, what) {
    pkg <- annotation(to)
    ok <- FALSE
    if (length(grep(".*db$", pkg))==1) {
        genes <- .toAnnotationDbi(tag, to, what)
        ok <- TRUE
    } else {
        if (!ok)
            tryCatch({
                genes <- .toAnnotationDirect(tag, to, what)
                ok <- TRUE
            }, error=function(err) {
                warning("direct map failed: ",  conditionMessage(err))
            })
        if (!ok)
            tryCatch({
                genes <- .toAnnotationRevmap(tag, to, what)
                ok <- TRUE
            }, error=function(err) {
                warning("reverse map failed: ", conditionMessage(err))
            })
    }
    if (!ok)
        .stopf("unable to map from '%s' to '%s'", tag, setType(to))
    genes
}

.toAnnotationDirect <- function(tag, to, what) {
    pkg <- annotation(to)
    map <- paste(pkg, tag, "2PROBE",sep="")
    .checkPackageAndMap(pkg, map)
    mapEnv <- get(map, envir=getNamespace(pkg))
    .getMappedGenes(genes(what), mapEnv, map, pkg)
}

.toAnnotationRevmap <- function(tag, to, what) {
    pkg <- annotation(to)
    map <- paste(gsub("db$", "", pkg), tag, sep="")
    .checkPackageAndMap(pkg, map)
    revMap <- revmap(get(map))
    genes <- genes(what)
    ogenes <- genes[sapply(genes, exists, envir=revMap)]
    .getMappedGenes(ogenes, revMap,
                    paste("revmap(", map, ")", sep=""), pkg)
}

.toAnnotationDbi <- function(tag, to, what) {
    pkg <- annotation(to)
    map <- paste(gsub("db$", "", pkg), tag, sep="")
    .checkPackageAndMap(pkg, map)
    mapEnv <- getNamespace(pkg)
    mapEnv <- revmap(get(map, envir=mapEnv))
    .getMappedGenes(genes(what), mapEnv, map, pkg)
}

setMethod("mapIdentifiers",
          ## this method resolves ambiguity between
          ## GeneSet#AnnotationIdentifier#GeneIdentifierType
          ## GeneSet#GeneIdentifierType#NullIdentifier
          ## for
          ## GeneSet#AnnotationIdentifier#NullIdentifier
          signature=signature(
            what="GeneSet",
            to="AnnotationIdentifier",
            from="NullIdentifier"),
          function(what, to, from, ...) {
              new(class(what), what, type=to)
          })

setMethod("mapIdentifiers",
          signature=signature(
            what="GeneSet",
            to="AnnotationIdentifier",
            from="GeneIdentifierType"), 
          function(what, from, to, ...) {
              tag <- toupper(setType(from))
              new(class(what), what,
                  genes=.toAnnotation(tag, to, what),
                  type=to)
          })
