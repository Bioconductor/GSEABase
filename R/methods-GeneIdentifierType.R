## GeneIdentifierType

.CONSTRUCTORS_GeneIdentifierType <- 
    names(getSubclasses(getClass("GeneIdentifierType")))

.CONSTRUCTORS_GeneIdentifierType <-
    .CONSTRUCTORS_GeneIdentifierType[!(.CONSTRUCTORS_GeneIdentifierType %in%
                                       c("AnnotationIdentifier"))]

.constructors_Simple(.CONSTRUCTORS_GeneIdentifierType[.CONSTRUCTORS_GeneIdentifierType!="AnnotationIdentifier"])

.constructors_Simple("AnnotationIdentifier", required="annotation")

.getters("GeneIdentifierType", c(geneIdType="type"))

setMethod("show",
          signature=signature(object="GeneIdentifierType"),
          function(object) cat("geneIdType:", geneIdType(object), "\n"))

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
              cat("geneIdType:", geneIdType(object),
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
          function(what, to, from, ..., verbose=FALSE) {
              callGeneric(what, to, from=geneIdType(what), ..., verbose=verbose)
          })

## Null --> X

setMethod("mapIdentifiers",
          signature=signature(
            what="GeneSet",
            to="GeneIdentifierType",
            from="NullIdentifier"),
          function(what, to, from, ..., verbose=FALSE) {
              new(class(what), what, geneIdType=to)
          })
          

## AnnotationIdentifier --> X

.checkPackageAndMap <- function(pkg, map) {
    if (!.requireQ(pkg))
        .stopf("cannot load annotation package '%s'", pkg)
    if (!exists(map, where=getNamespace(pkg)))
        .stopf("map '%s' not found in annotation package '%s'",
               map, pkg)
}

.getMappedGenes <- function(geneIds, mapEnv, map, pkg, verbose=FALSE) {
    ngenes <- mget(geneIds, mapEnv, ifnotfound=as.character(NA))
    if (verbose && (any(length(ngenes) != 1) || any(is.na(ngenes))))
        .warningf("annotation map '%s' is not 1:1 in '%s'",
                  map, pkg)
    ugenes <- unique(unlist(ngenes))
    if (is.null(ugenes))
        character(0)
    else
        as.character(ugenes[!is.na(ugenes)])
}

## tag: e.g., ENTREZID
.fromAnnotation <- function(from, tag, what, verbose=FALSE) {
    pkg <- annotation(from)
    map <- paste(gsub(".db$", "", pkg), tag, sep="")
    .checkPackageAndMap(pkg, map)
    mapEnv <-
        get(map, envir=as.environment(paste("package", pkg, sep=":")),
            inherits=FALSE)
    .getMappedGenes(geneIds(what), mapEnv, map, pkg, verbose=verbose)
}

setMethod("mapIdentifiers",
          signature=signature(
            what="GeneSet",
            to="GeneIdentifierType",
            from="AnnotationIdentifier"),
          function(what, to, from, ..., verbose=FALSE) {
              tag <- toupper(geneIdType(to))
              new(class(what), what,
                  geneIds=.fromAnnotation(from, tag, what, verbose=verbose),
                  geneIdType=to)
          })

## X --> AnnotationIdentifier

## tag: e.g., SYMBOL
.toAnnotation <- function(tag, to, what, verbose=FALSE) {
    pkg <- annotation(to)
    ok <- FALSE
    if (length(grep(".*db$", pkg))==1) {
        geneIds <- .toAnnotationDbi(tag, to, what, verbose=verbose)
        ok <- TRUE
    } else {
        if (!ok)
            tryCatch({
                geneIds <- .toAnnotationDirect(tag, to, what, verbose=verbose)
                ok <- TRUE
            }, error=function(err) {
                if (verbose)
                    warning("direct map failed: ",  conditionMessage(err))
            })
        if (!ok)
            tryCatch({
                geneIds <- .toAnnotationRevmap(tag, to, what, verbose=verbose)
                ok <- TRUE
            }, error=function(err) {
                if (verbose)
                    warning("reverse map failed: ", conditionMessage(err))
            })
    }
    if (!ok)
        .stopf("unable to map from '%s' to '%s'", tag, geneIdType(to))
    geneIds
}

.toAnnotationDirect <- function(tag, to, what, verbose=FALSE) {
    pkg <- annotation(to)
    map <- paste(pkg, tag, "2PROBE",sep="")
    .checkPackageAndMap(pkg, map)
    mapEnv <-
        get(map,
            envir=as.environment(paste("package", pkg, sep=":")),
            inherits=FALSE)
    .getMappedGenes(geneIds(what), mapEnv, map, pkg, verbose=verbose)
}

.toAnnotationRevmap <- function(tag, to, what, verbose=FALSE) {
    pkg <- annotation(to)
    map <- paste(gsub(".db$", "", pkg), tag, sep="")
    .checkPackageAndMap(pkg, map)
    revMap <-
        revmap(get(map,
                   envir=as.environment(paste("package", pkg, sep=":")),
                   inherits=FALSE))
    geneIds <- geneIds(what)
    ogenes <- geneIds[sapply(geneIds, exists, envir=revMap)]
    .getMappedGenes(ogenes, revMap,
                    paste("revmap(", map, ")", sep=""), pkg, verbose=verbose)
}

.toAnnotationDbi <- function(tag, to, what, verbose=FALSE) {
    pkg <- annotation(to)
    map <- paste(gsub(".db$", "", pkg), tag, sep="")
    .checkPackageAndMap(pkg, map)
    mapEnv <- getNamespace(pkg)
    mapEnv <- revmap(get(map, envir=mapEnv, inherits=FALSE))
    .getMappedGenes(geneIds(what), mapEnv, map, pkg, verbose=verbose)
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
          function(what, to, from, ..., verbose=FALSE) {
              new(class(what), what, geneIdType=to)
          })

setMethod("mapIdentifiers",
          signature=signature(
            what="GeneSet",
            to="AnnotationIdentifier",
            from="GeneIdentifierType"), 
          function(what, from, to, ..., verbose=FALSE) {
              tag <- toupper(geneIdType(from))
              new(class(what), what,
                  geneIds=.toAnnotation(tag, to, what, verbose=verbose),
                  geneIdType=to)
          })
