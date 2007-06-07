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

## construct these programatically?

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
        stop(sprintf("cannot load annotation package '%s'", pkg))
    if (!exists(map, where=getNamespace(pkg)))
        stop(sprintf("map '%s' not found in annotation package '%s'",
                     map, pkg))
}

.getMappedGenes <- function(genes, map, mapEnv, pkg, getter=mget) {
    ngenes <- getter(genes, mapEnv)
    if (any(length(ngenes) != 1) || any(is.na(ngenes)))
        warning(sprintf("annotation map '%s' is not 1:1 in '%s'",
                        map, pkg))
    ugenes <- unique(unlist(ngenes))
    as.character(ugenes[!is.na(ugenes)])
}

## tag: e.g., ENTREZID
.fromAnnotation <- function(from, tag, what) {
    pkg <- annotation(from)
    map <- paste(gsub("db$", "", pkg), tag, sep="")
    .checkPackageAndMap(pkg, map)
    mapEnv <- get(map, envir=getNamespace(pkg))
    .getMappedGenes(genes(what), map, mapEnv, pkg)
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
        stop(sprintf("unable to map from '%s' to '%s'",
                     tag, setType(to)))
    genes
}

.toAnnotationDirect <- function(tag, to, what) {
    pkg <- annotation(to)
    map <- paste(pkg, tag, "2PROBE",sep="")
    .checkPackageAndMap(pkg, map)
    mapEnv <- get(map, envir=getNamespace(pkg))
    .getMappedGenes(genes(what), map, mapEnv, pkg)
}

.toAnnotationRevmap <- function(tag, to, what) {
    pkg <- annotation(to)
    map <- paste(gsub("db$", "", pkg), tag, sep="")
    .checkPackageAndMap(pkg, map)
    mapEnv <- get(map)
    wgenes <- genes(what)
    revMap <- new.env(hash=TRUE, parent=emptyenv(),
                      size=length(genes))
    probeIds <- ls(mapEnv)
    entrezIds <- mget(probeIds, mapEnv)
    okIds <- lapply(entrezIds, match, wgenes, nomatch=0)
    okIds <- lapply(okIds, `>`, 0)
    mapply(function(probeId, entrezIds, okIds) {
        if (any(okIds)) {
            entrezIds <- entrezIds[okIds]
            for (eid in entrezIds)
                if (exists(eid, revMap))
                    revMap[[eid]] <- c(revMap[[eid]], probeId)
                else
                    revMap[[eid]] <- probeId
        }
    }, probeIds, entrezIds, okIds)
    if (length(revMap) != length(wgenes) ||
        any(eapply(revMap, length) != 1))
        warning(sprintf("annotation map '%s' is not 1:1 in '%s'",
                        map, pkg))
    as.character(unique(unlist(mget(ls(revMap), revMap))))
}

.toAnnotationDbi <- function(tag, to, what) {
    pkg <- annotation(to)
    map <- paste(gsub("db$", "", pkg), tag, sep="")
    .checkPackageAndMap(pkg, map)
    mapEnv <- getNamespace(pkg)
    mapEnv <- revmap(get(map, envir=mapEnv))
    .getMappedGenes(genes(what), map, mapEnv, pkg, getter=AnnotationDbi:::mget)
}

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
