## Constructors

setMethod("GeneSetCollection",
          signature=signature(
            object="GeneSet",
            idType="missing",
            setType="missing"),
          function(object, ..., idType, setType) {
              new("GeneSetCollection", list(object, ...))
          })

setMethod("GeneSetCollection",
          signature=signature(
            object="list",
            idType="missing",
            setType="missing"),
          function(object, ..., idType, setType) {
              names(object) <- NULL
              new("GeneSetCollection", object)
          })

## Use revmap

.GSC_filter_by_probe <- function(genes, probes) {
    probesOk <- lapply(genes, "%in%", probes)
    genes <- mapply("[", genes, probesOk)
    genes[sapply(genes, length) != 0]
}

.GSC_genes_helper <- function(idType, setType) {
    mapEnv <- revmap(getAnnMap(toupper(collectionType(setType)),
                               annotation(idType)))
    lapply(as.list(mapEnv), unique)
}

.GSC_ExpressionSet_helper <- function(object, idType, setType) {
    .GSC_filter_by_probe(.GSC_genes_helper(idType, setType),
                        featureNames(object))
}

.GSC_CollectionIdTypes <- function(genes, setType) {
    collType <- paste(collectionType(setType), "Collection",
                      sep="")
    lapply(names(genes),
           function(ids) do.call(collType, list(ids=ids)))
}

.GSC_CollectionType <- function(genes, idType, collTypes) {
    organism <- getAnnMap("ORGANISM", annotation(idType))
    gss <- mapply(GeneSet,
                  genes,
                  setName=names(genes),
                  collectionType=collTypes,
                  MoreArgs=list(
                    geneIdType=idType,
                    organism=organism))
    GeneSetCollection(gss)
}

setMethod("GeneSetCollection",
          signature=signature(
            object="character",
            idType="AnnotationIdentifier",
            setType="CollectionType"),
          function(object, ..., idType, setType) {
              helper <- .GSC_genes_helper(idType, setType)
              genes <- .GSC_filter_by_probe(helper, object)
              collTypes <- .GSC_CollectionIdTypes(genes, setType)
              .GSC_CollectionType(genes, idType, collTypes)
          })

setMethod("GeneSetCollection",
          signature=signature(
            object="character",
            idType="AnnotationIdentifier",
            setType="CollectionIdType"),
          function(object, ..., idType, setType) {
              helper <- .GSC_genes_helper(idType, setType)
              genes <- .GSC_filter_by_probe(helper, object)
              if (length(ids(setType)) > 0)
                  genes <- genes[names(genes) %in% ids(setType)]
              collTypes <- .GSC_CollectionIdTypes(genes, setType)
              .GSC_CollectionType(genes, idType, collTypes)
          })

setMethod("GeneSetCollection",
          signature=signature(
            object="ExpressionSet",
            idType="missing",
            setType="CollectionType"),
          function(object, ..., idType, setType) {
              idType <- AnnotationIdentifier(annotation(object))
              genes <- .GSC_ExpressionSet_helper(object, idType, setType)
              collTypes <- lapply(names(genes), setType)
              .GSC_CollectionType(genes, idType, collTypes)
          })

setMethod("GeneSetCollection",
          signature=signature(
            object="ExpressionSet",
            idType="missing",
            setType="CollectionIdType"),
          function(object, ..., idType, setType) {
              idType <- AnnotationIdentifier(annotation(object))
              genes <- .GSC_ExpressionSet_helper(object, idType, setType)
              if (length(ids(setType)) > 0)
                  genes <- genes[names(genes) %in% ids(setType)]
              collTypes <- .GSC_CollectionIdTypes(genes, setType)
              .GSC_CollectionType(genes, idType, collTypes)
          })

setMethod("GeneSetCollection",
          signature=signature(
            object="missing",
            idType="AnnotationIdentifier",
            setType="CollectionType"),
          function(object, ..., idType, setType) {
              genes <- .GSC_genes_helper(idType, setType)
              .GSC_CollectionType(genes, idType, setType)
          })

setMethod("GeneSetCollection",
          signature=signature(
            object="missing",
            idType="AnnotationIdentifier",
            setType="CollectionIdType"),
          function(object, ..., idType, setType) {
              genes <- .GSC_genes_helper(idType, setType)
              if (length(ids(setType))>0)
                  genes <- genes[names(genes) %in% ids(setType)]
              collTypes <- .GSC_CollectionIdTypes(genes, setType)
              .GSC_CollectionType(genes, idType, collTypes)
          })

## Use direct map

setMethod("GeneSetCollection",
          signature=signature(
            object="ExpressionSet",
            idType="missing",
            setType="KEGGCollection"),
          function(object, ..., idType, setType) {
              idType <- AnnotationIdentifier(annotation(object))
              genes <- as.list(getAnnMap("PATH2PROBE", annotation(idType)))
              genes <- .GSC_filter_by_probe(genes, featureNames(object))
              if (length(ids(setType))>0)
                  genes <- genes[names(genes) %in% ids(setType)]
              collTypes <- .GSC_CollectionIdTypes(genes, setType)
              .GSC_CollectionType(genes, idType, collTypes, ...)
          })

setMethod("GeneSetCollection",
          signature=signature(
            object="missing",
            idType="AnnotationIdentifier",
            setType="KEGGCollection"),
          function(object, ..., idType, setType) {
              genes <- as.list(getAnnMap("PATH2PROBE", annotation(idType)))
              if (length(ids(setType)) > 0)
                  genes <- genes[names(genes) %in% ids(setType)]
              collTypes <- .GSC_CollectionIdTypes(genes, setType)
              .GSC_CollectionType(genes, idType, collTypes, ...)
          })

.GSC_GO_helper <- function(genes, idType, setType, ...) {
    ## filter on evidence codes
    evidenceCode <- evidenceCode(setType)
    eviOk <- lapply(lapply(genes, names), "%in%", evidenceCode)
    genes <- mapply("[", genes, eviOk)
    ugenes <- lapply(genes, unique)
    ugenes <- ugenes[sapply(ugenes, length) != 0]

    organism <- getAnnMap("ORGANISM", annotation(idType))
    gss <- mapply(function(ids, setName, collectionType, ...) {
        GeneSet(ids,
                setName=setName,
                collectionType=GOCollection(
                  ids=setName,
                  evidenceCode=evidenceCode(collectionType),
                  ontology=ontology(collectionType)),
                ...)
    }, ugenes, setName=names(ugenes), MoreArgs=list(
                                collectionType=setType,
                                geneIdType=idType,
                                organism=organism, ...))
    GeneSetCollection(gss)
}

setMethod("GeneSetCollection",
          signature=signature(
            object="ExpressionSet",
            idType="missing",
            setType="GOCollection"),
          function(object, ..., idType, setType) {
              idType <- AnnotationIdentifier(annotation(object))
              genes <- as.list(getAnnMap("GO2PROBE", annotation(idType)))
              genes <- .GSC_filter_by_probe(genes, featureNames(object))
              .GSC_GO_helper(genes, idType=idType, setType=setType, ...)
          })

setMethod("GeneSetCollection",
          signature=signature(
            object="missing",
            idType="AnnotationIdentifier",
            setType="GOCollection"),
          function(object, ..., idType, setType) {
              genes <- as.list(getAnnMap("GO2PROBE", annotation(idType)))
              .GSC_GO_helper(genes, idType=idType, setType=setType, ...)
          })

## updateObject

setMethod("updateObject",
          signature=signature(
            object="GeneSetCollection"),
          function(object, ..., verbose=FALSE) {
              if (verbose)
                  message("updateObject,GeneSetCollection-method")
              initialize(object,
                         lapply(object, updateObject, verbose=verbose))
          })

## accessors

setMethod("geneIds<-",
          signature=signature(
            object="GeneSetCollection",
            value="list"),
          function(object, value) {
              lapply(object, "geneIds<-", value)
              object
          })

setMethod("geneIds",
          signature=signature(
            object="GeneSetCollection"),
          function(object) {
              lapply(object, geneIds)
          })

setMethod("names",
          signature=signature(
            x="GeneSetCollection"),
          function(x) sapply(x, setName))

## [, [[

.characterToIndex <- function(x, i) {
    if (any(duplicated(i)))
        .stopf("duplicate setNames not allowed: '%s'",
               paste(i[duplicated(i)], collapse="', '"))
    idx <- pmatch(i, names(x))
    if (any(is.na(idx)))
        .stopf("unknown setNames: '%s'",
               paste(i[is.na(idx)], collapse="', '"))
    idx
}

.subset <- function(x, i) {
    x@.Data <- x@.Data[i]
    x
}

setMethod("[",
          signature=signature(
            x="GeneSetCollection",
            i="logical"),
          function(x, i, j, ..., drop=TRUE) {
              if (length(i) > length(x))
                  .stopf("logical length '%d' greater than GeneSetCollection length '%d'",
                         length(i), length(x))
              .subset(x, i)
          })

setMethod("[",
          signature=signature(
            x="GeneSetCollection",
            i="numeric"),
          function(x, i, j, ..., drop=TRUE) {
              if (any(duplicated(i)))
                  .stopf("duplicate index not allowed: '%s'",
                         paste(i[duplicated(i)], collapse="', '"))
              if (any(i > length(x)))
                  .stopf("subscript out of bounds: '%s'",
                         paste(i[i>length(x)], collapse="', '"))
              .subset(x, i)
          })

setMethod("[",
          signature=signature(
            x="GeneSetCollection",
            i="character"),
          function(x, i, j, ..., drop=TRUE) {
              idx <- .characterToIndex(x, i)
              .subset(x, idx)
          })

setMethod("[[",
          signature=signature(
            x="GeneSetCollection",
            i="character"),
          function(x, i, j, ...) {
              idx <- match(i, names(x))
              if (length(i) != 1 || is.na(idx))
                  .stopf("subscript out of bounds: '%s'",
                         paste(i, collapse="', "))
              x[[idx]]
          })

## [<-, [[<-

.subsetReplace <- function(x, i, value) {
    x@.Data[i] <- value
    x
}

setReplaceMethod("[",
                 signature=signature(
                   x="GeneSetCollection",
                   value="ANY"),
                 function(x, i, j, ..., value) {
                     .stopf("cannot assign object of class '%s' to '%s'",
                            class(value), class(x))
                 })

setReplaceMethod("[",
                 signature=signature(
                   x="GeneSetCollection",
                   value="GeneSet"),
                 function(x, i, j, ..., value) {
                     .subsetReplace(x, i, value)
                 })

setReplaceMethod("[",
                 signature=signature(
                   x="GeneSetCollection",
                   i="character",
                   value="GeneSet"),
                 function(x, i, j, ..., value) {
                     idx <- .characterToIndex(x, i)
                     .subsetReplace(x, idx, value)
                 })

setReplaceMethod("[[",
                 signature=signature(
                   x="GeneSetCollection",
                   value="ANY"),
                 function(x, i, j ,..., value) {
                     .stopf("cannot assign object of class '%s' to '%s'",
                            class(value), class(x))
                 })

setReplaceMethod("[[",
                 signature=signature(
                   x="GeneSetCollection",
                   i="numeric",
                   value="GeneSet"),
                 function(x, i, j ,..., value) {
                     .subsetReplace(x, i, value)
                 })

setReplaceMethod("[[",
                 signature=signature(
                   x="GeneSetCollection",
                   i="character",
                   value="GeneSet"),
                 function(x, i, j, ..., value) {
                     if (length(i)!=1)
                         .stopf("index must be length 1, but is '%d'",
                                length(i))
                     idx <- match(i, names(x))
                     if (is.na(idx))
                         .stopf("only replacement of existing setNames supported")
                     .subsetReplace(x, idx, value)
                 })

## logic

setMethod("&",
          signature=signature(
            e1="GeneSetCollection",
            e2="ANY"),
          function(e1, e2) {
              GeneSetCollection(lapply(e1, `&`, e2))
          })

setMethod("&",
          signature=signature(
            e1="ANY",
            e2="GeneSetCollection"),
          ## FIXME: Should be able to define the complement on Logic,
          ## but this does not work 2007-08-14 r42505
          function(e1, e2) e2 & e1)

setMethod("&",
          signature=signature(
            e1="GeneSetCollection",
            e2="GeneSetCollection"),
          function(e1, e2) {
              .stopf("'%s' & '%s' not yet implemented",
                     class(e1), class(e2))
          })

setMethod("intersect",
          signature=signature(
            x="GeneSetCollection",
            y="ANY"),
          function(x, y) x & y)

setMethod("intersect",
          signature=signature(
            x="ANY",
            y="GeneSetCollection"),
          function(x, y) y & x)

setMethod("|",
          signature=signature(
            e1="GeneSetCollection",
            e2="ANY"),
          function(e1, e2) {
              GeneSetCollection(lapply(e1, `|`, e2))
          })

setMethod("|",
          signature=signature(
            e1="ANY",
            e2="GeneSetCollection"),
          ## FIXME: Should be able to define the complement on Logic,
          ## but this does not work 2007-08-14 r42505
          function(e1, e2) e2 | e1)

setMethod("|",
          signature=signature(
            e1="GeneSetCollection",
            e2="GeneSetCollection"),
          function(e1, e2) {
              .stopf("'%s' | '%s' not yet implemented",
                     class(e1), class(e2))
          })

setMethod("union",
          signature=signature(
            x="GeneSetCollection",
            y="ANY"),
          function(x, y) x | y)

setMethod("union",
          signature=signature(
            x="ANY",
            y="GeneSetCollection"),
          function(x, y) y | x)

setMethod("setdiff",
          signature=signature(
            x="GeneSetCollection",
            y="ANY"),
          function(x, y) {
              GeneSetCollection(lapply(x, setdiff, y))
          })

## setMethod("setdiff",
##           signature=signature(
##             x="GeneSetCollection",
##             y="GeneSetCollection"),
##           function(x, y) {
##               .stopf("'setdiff(%s, %s)' not yet implemented",
##                      class(x), class(y))
##           })

setMethod("Logic",
          signature=signature(e1="character", e2="GeneSetCollection"),
          function(e1, e2) callGeneric(e2, e1))

setMethod("Logic",
          signature=signature(e1="GeneSet", e2="GeneSetCollection"),
          function(e1, e2) callGeneric(e2, e1))

## mapIdentifiers

setMethod("mapIdentifiers",
          signature=signature(
            what="GeneSetCollection",
            to="GeneIdentifierType",
            from="missing"),
          function(what, to, from, ..., verbose=FALSE) {
              GeneSetCollection(lapply(what, mapIdentifiers, to, ...,
                                       verbose=verbose))
          })

## incidence

setMethod("incidence",
          signature=signature(
            x="GeneSetCollection"),
          function(x, ...) {
              args <- c(x, ...)
              .incidence(lapply(args, geneIds),
                         lapply(args, setName))
          })

## show

setMethod("show",
          signature=signature(
            object="GeneSetCollection"),
          function(object) {
              some <- function(x)
                  paste(paste(Biobase::selectSome(x, 4), collapse=", "),
                        " (", length(x), " total)", sep="")
              gids <- unique(unlist(geneIds(object)))
              itypes <- unique(sapply(lapply(object, geneIdType), class))
              ctypes <- unique(sapply(lapply(object, collectionType), class))
              cat("GeneSetCollection\n",
                  "  names: ", some(names(object)), "\n",
                  "  unique identifiers: ", some(gids), "\n",
                  "  types in collection:\n",
                  "    geneIdType: ", some(itypes), "\n",
                  "    collectionType: ", some(ctypes), "\n",
                  sep="")
          })
