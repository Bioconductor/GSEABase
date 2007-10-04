## constructors / getters / setters in AllClasses.R

IdFactory <- function(classPrefix,
                      type=classPrefix,
                      where=topenv(parent.frame()), ...) {
    class <- paste(classPrefix, "Identifier", sep="")
    eval(substitute({
        setClass(CLASS,
                 contains = "GeneIdentifierType",
                 prototype = prototype(type=TYPE),
                 where=WHERE)
        .constructors_Simple(CLASS, where=WHERE)
        if (existsMethod("initialize",
                         signature=CLASS, where=WHERE))
            removeMethod("initialize",
                         signature=CLASS, where=WHERE)
        if (existsMethod("show",
                         signature=CLASS, where=WHERE))
            removeMethod("show",
                         signature=CLASS, where=WHERE)
    }, list(CLASS=class,
            TYPE=mkScalar(type),
            WHERE=where)))
    invisible(get(class, envir=where))
}

## GeneIdentifierType

setMethod("show",
          signature=signature(object="GeneIdentifierType"),
          function(object) cat("geneIdType:", geneIdType(object), "\n"))

## AnnotationIdentifier

AnnotationIdFactory <- function(classPrefix,
                                type=classPrefix,
                                where=topenv(parent.frame()), ...) {
    class <- paste(classPrefix, "Identifier", sep="")
    eval(substitute({
        setClass(CLASS,
                 contains = "AnnotationIdentifier",
                 prototype = prototype(type=TYPE),
                 where = WHERE)
        GSEABase:::.constructors_Simple(CLASS,
                                        required="annotation", where=WHERE)
        setMethod("initialize",
                  signature=signature(.Object=CLASS),
                  function(.Object, .Template=.Object, ...,
                           annotation = .Template@annotation) {
                      callNextMethod(.Object, .Template, ...,
                                     annotation=mkScalar(annotation))
                  }, where=WHERE)
        setMethod("show",
                  signature=signature(object=CLASS),
                  function(object) {
                      cat("geneIdType:", geneIdType(object),
                          if (nzchar(annotation(object))) {
                              paste("(", annotation(object), ")", sep="")
                          }, "\n")
                  }, where=WHERE)
    }, list(CLASS=class,
            TYPE=mkScalar(type),
            WHERE=where)))
    invisible(get(class, envir=where))
}

setMethod("initialize",
          signature=signature(.Object="AnnotationIdentifier"),
          function(.Object, .Template=.Object, ...,
                   annotation=.Template@annotation) {
              callNextMethod(.Object, .Template, ...,
                             annotation = mkScalar(annotation))
          })

.getAnnMap <- function(object, symbol) {
    pkgName <- annotation(object)
    pkg <- sub(".db$", "", pkgName)
    symbol <- paste(pkg, symbol, sep="")
    if (!.requireQ(pkgName)) {
        pkgNameDb <- paste(pkg, "db", sep=".")
        if (pkgName != pkgNameDb && !.requireQ(pkgNameDb))
            .stopf("cannot load annotation package '%s'", pkgName)
        pkgName <- pkgNameDb
    }
    pkgPos <- match(paste("package", pkgName, sep=":"), search())
    if (!exists(symbol, pkgPos, inherits=FALSE))
        .stopf("no symbol '%s' in annotation package '%s'",
               symbol, pkgName)
    get(symbol, pkgPos, inherits=FALSE)
}

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
              callGeneric(what, to, from=geneIdType(what), ...,
                          verbose=verbose)
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

.getMappedGenes <- function(geneIds, mapEnv, map, pkg, verbose=FALSE) {
    anNA <- function(x) length(x) ==1 && is.na(x)
    ngenes <- mget(geneIds, mapEnv, ifnotfound=as.character(NA))
    if (verbose) {
        if (any(sapply(ngenes, length) != 1))
            .warningf("annotation map '%s' not 1:1 in '%s'\n  ids: '%s'",
                      map, pkg,
                      paste(names(ngenes)[sapply(ngenes, length) != 1],
                            collapse="', '"))
        if (any(sapply(ngenes, anNA)))
            .warningf("annotation map '%s' had %d 'NA' values in '%s'",
                      map, sum(sapply(ngenes, anNA)), pkg)
    }
    ugenes <- unique(unlist(ngenes))
    if (verbose && (length(ugenes) != length(geneIds)))
        .warningf("annotation map '%s' is %d:%d (not 1:1) in '%s'",
                  map, length(geneIds), length(ugenes), pkg)
    if (is.null(ugenes))
        character(0)
    else
        as.character(ugenes[!is.na(ugenes)])
}

## tag: e.g., ENTREZID
.fromAnnotation <- function(from, tag, what, verbose=FALSE) {
    pkg <- annotation(from)
    map <- paste(gsub(".db$", "", pkg), tag, sep="")
    mapEnv <- .getAnnMap(from, tag)
    .getMappedGenes(geneIds(what), mapEnv, map, pkg, verbose=verbose)
}

setMethod("mapIdentifiers",
          signature=signature(
            what="GeneSet",
            to="GeneIdentifierType",
            from="AnnotationIdentifier"),
          function(what, to, from, ..., verbose=FALSE) {
              tag <- toupper(geneIdType(to))
              geneIds <- .fromAnnotation(from, tag, what,
                                         verbose=verbose)
              new(class(what), what,
                  geneIds=geneIds, geneIdType=to)
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
    mapEnv <- .getAnnMap(to, paste(tag, "2PROBE", sep=""))
    .getMappedGenes(geneIds(what), mapEnv, map, pkg, verbose=verbose)
}

.toAnnotationRevmap <- function(tag, to, what, verbose=FALSE) {
    pkg <- annotation(to)
    map <- paste(gsub(".db$", "", pkg), tag, sep="")
    revMap <- revmap(.getAnnMap(to, tag))
    geneIds <- geneIds(what)
    ogenes <- geneIds[sapply(geneIds, exists, envir=revMap)]
    .getMappedGenes(ogenes, revMap,
                    paste("revmap(", map, ")", sep=""), pkg, verbose=verbose)
}

.toAnnotationDbi <- function(tag, to, what, verbose=FALSE) {
    pkg <- annotation(to)
    map <- paste(gsub(".db$", "", pkg), tag, sep="")
    mapEnv <- revmap(.getAnnMap(to, tag))
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
          function(what, to, from, ..., verbose=FALSE) {
              tag <- toupper(geneIdType(from))
              new(class(what), what,
                  geneIds=.toAnnotation(tag, to, what, verbose=verbose),
                  geneIdType=to)
          })

setMethod("mapIdentifiers",
          signature=signature(
            what="GeneSet",
            to="GeneIdentifierType",
            from="environment"),
          function(what, to, from, ..., verbose=FALSE) {
              geneIds <- .getMappedGenes(geneIds(what),
                                         from,
                                         "environment",
                                         "user-supplied environment",
                                         verbose=verbose)
              new(class(what), what,
                  geneIds = geneIds,
                  geneIdType = to)
          })
