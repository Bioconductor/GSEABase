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

.GeneIdType_asString <- function(x) {
    paste(geneIdType(x),
          if (nchar(annotation(x)) > 0) {
              paste(" (", annotation(x), ")", sep="")
          }, sep="")
}

setMethod("show",
          signature=signature(object="GeneIdentifierType"),
          function(object) {
              cat("geneIdType: ", .GeneIdType_asString(object), "\n",
                  sep="")
          })

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
                      cat("geneIdType: ", .GeneIdType_asString(object),
                          "\n", sep="")
                  }, where=WHERE)
    }, list(CLASS=class,
            TYPE=mkScalar(type),
            WHERE=where)))
    invisible(get(class, envir=where))
}

setMethod("initialize",
          signature=signature(.Object="GeneIdentifierType"),
          function(.Object, .Template=.Object, ...,
                   annotation=.Template@annotation) {
              callNextMethod(.Object, .Template, ...,
                             annotation = mkScalar(annotation))
          })

.getAnnMap <- function(object, symbol) {
    getAnnMap(symbol, annotation(object))
}

setMethod("show",
          signature=signature(object="AnnotationIdentifier"),
          function(object) {
              cat("geneIdType: ", .GeneIdType_asString(object), "\n",
                  sep="")
          })

## mapIdentifiers

.mapString <- function(from, to) {
    paste(.GeneIdType_asString(from), "->",
          .GeneIdType_asString(to))
}
    

.mapIdentifiers_isNullMap <- function(from, to, verbose) {
    res <- (nchar(annotation(from)) > 0 &&
            nchar(annotation(to)) > 0 &&
            geneIdType(from) == geneIdType(to))
    if (verbose && res) {
        .warningf("map '%s' is an identity; no map to perform",
                  .mapString(from, to))
    }
    res
}

.mapIdentifiers_isMappable <- function(from, to) {
    isDifferentAnno <- (nchar(annotation(from)) > 0 &&
                        nchar(annotation(to)) > 0 &&
                        annotation(from) != annotation(to))
    if (isDifferentAnno)
      reason <- "'from' and 'to' annotations differ"
    isBothUnnamed <- (!isDifferentAnno &&
                      nchar(annotation(from)) == 0 &&
                      nchar(annotation(to)) == 0)
    if (isBothUnnamed)
      reason <- "neither GeneIdentifierType has annotation"
    if (isDifferentAnno || isBothUnnamed)
      .stopf("unable to map from '%s' to '%s'\n    %s",
             .GeneIdType_asString(from),
             .GeneIdType_asString(to), reason)
    TRUE
}

.mapIdentifiers_normalize <- function(from, to) {
    .mapIdentifiers_isMappable(from, to)
    if (nchar(annotation(from)) == 0)
      annotation(from) <- annotation(to)
    else if (nchar(annotation(to)) == 0)
      annotation(to) <- annotation(from)
    list(from, to)
}


.mapIdentifiers_selectMaps <- function(from, to) {
    isOrg <- function(x) {
        ## is this an 'org' package??
        length(grep("^org\\.", annotation(x)))==1
    }
    if (is(from, "AnnotationIdentifier")) {
        ## one map: AnnotationIdentifier --> to
        first <- getAnnMap(toupper(geneIdType(to)), annotation(from))
        second <- NULL
    } else if (is(to, "AnnotationIdentifier")) {
        ## one map: revmap(AnnotationIdentifier --> from)
        map <- getAnnMap(toupper(geneIdType(from)), annotation(to))
        first <- revmap(map)
        second <- NULL
    } else if (is(from, "EntrezIdentifier") && isOrg(from)) {
        ## one map: EntrezIdentifier --> to
        first <- getAnnMap(toupper(geneIdType(to)), annotation(from))
        second <- NULL
    } else if (is(to, "EntrezIdentifier") && isOrg(to)) {
        ## one map: revmap(EntrezIdentifier --> to)
        map <- getAnnMap(toupper(geneIdType(from)), annotation(to))
        first <- revmap(map)
        second <- NULL
    } else {
        ## two maps
        map <- getAnnMap(toupper(geneIdType(from)), annotation(from))
        first <- revmap(map)
        second <- getAnnMap(toupper(geneIdType(to)), annotation(to))
    }
    c(first=first, second=second)
}

.mapIdentifiers_doWithMap <- function(keys, map, from, to, verbose) {
    isNA <- function(x) length(x)==1 && is.na(x)
    vals <- mget(keys, map, ifnotfound=NA_character_)
    if (verbose) {
        if (any(sapply(vals, length) != 1))
          .warningf("map '%s' not 1:1\n  ids: '%s'",
                    .mapString(from, to),
                    paste(names(vals)[sapply(vals, length) != 1],
                          collapse="', '"))
        if (any(sapply(vals, isNA)))
          .warningf("map '%s' had %d 'NA' values",
                    .mapString(from, to), sum(sapply(vals, isNA)))
    }
    uvals <- unique(unlist(vals))
    if (verbose && (length(uvals) != length(keys)))
      .warningf("map '%s' is %d:%d (not 1:1)",
                .mapString(from, to), length(keys), length(uvals))
    if (is.null(uvals))
      character(0)
    else
      as.character(uvals[!is.na(uvals)])
}

.mapIdentifiers_map <- function(ids, from, to, verbose=FALSE) {
    doMap <- .mapIdentifiers_doWithMap  # abbrevation; nothing fancy
    map <- .mapIdentifiers_selectMaps(from, to)
    if (length(map)==1) doMap(ids, map[[1]], from, to, verbose)
    else {
        key <- doMap(ids, map[[1]], from, to, verbose)
        doMap(key, map[[2]], from, to, verbose)
    }
}

## Methods

setMethod("mapIdentifiers",
          signature=signature(
            what="GeneSet",
            to="GeneIdentifierType",
            from="missing"),
          function(what, to, from, ..., verbose=FALSE) {
              callGeneric(what, to, from=geneIdType(what), ...,
                          verbose=verbose)
          })

setMethod("mapIdentifiers",
          signature=signature(
            what="GeneSet",
            to="GeneIdentifierType",
            from="NullIdentifier"),
          function(what, to, from, ..., verbose=FALSE) {
              ## always valid
              initialize(what, geneIdType=to)
          })

setMethod("mapIdentifiers",
          signature=signature(
            what="GeneSet",
            to="NullIdentifier",
            from="GeneIdentifierType"),
          function(what, to, from, ..., verbose=FALSE) {
              ## always valid
              initialize(what, geneIdType=to)
          })
          
setMethod("mapIdentifiers",
          signature=signature(
            what="GeneSet",
            to="GeneIdentifierType",
            from="GeneIdentifierType"),
          function(what, to, from, ..., verbose=FALSE) {
              type <- .mapIdentifiers_normalize(from, to)
              if (.mapIdentifiers_isNullMap(type[[1]], type[[2]],
                                            verbose))
                  return(what)
              ids <- geneIds(what)
              ids <- .mapIdentifiers_map(ids, type[[1]], type[[2]],
                                         verbose)
              initialize(what, geneIds=ids, geneIdType=type[[2]])
          })

setMethod("mapIdentifiers",
          signature=signature(
            what="GeneSet",
            to="GeneIdentifierType",
            from="environment"),
          function(what, to, from, ..., verbose=FALSE) {
              doMap <- .mapIdentifiers_doWithMap
              ids <- doMap(geneIds(what), from,
                           "environment", "user-supplied environment",
                           verbose=verbose)
              initialize(what, geneIds=ids, geneIdType=to)
          })

setMethod("mapIdentifiers",
          signature=signature(
            what="GeneSet",
            to="GeneIdentifierType",
            from="AnnDbBimap"),
          function(what, to, from, ..., verbose=FALSE) {
              doMap <- .mapIdentifiers_doWithMap
              ids <- doMap(geneIds(what), from,
                           deparse(substitute(from)),
                           "user-supplied AnnDbBimap", verbose=verbose)
              initialize(what, geneIds=ids, geneIdType=to)
          })
