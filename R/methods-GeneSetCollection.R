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
              new("GeneSetCollection", object)
          })

.GSC_KEGG_helper <- function(genes, idType, setType, ...)
{
    organism <- .getAnnMap(idType, "ORGANISM")
    gss <- mapply(function(genes, setName, ...) {
        GeneSet(genes, setName=setName, ...)
    }, genes, names(genes), MoreArgs=list(
                              collectionType=setType,
                              geneIdType=idType,
                              organism=organism))

    GeneSetCollection(gss)
}

.GSC_filter_by_probe <- function(genes, probes) {
    probesOk <- lapply(genes, "%in%", probes)
    genes <- mapply("[", genes, probesOk)
    genes[sapply(genes, length) != 0]
}

setMethod("GeneSetCollection",
          signature=signature(
            object="missing",
            idType="AnnotationIdentifier",
            setType="KEGGCollection"),
          function(object, ..., idType, setType) {
              genes <- as.list(.getAnnMap(idType, "PATH2PROBE"))
              .GSC_KEGG_helper(genes, idType, setType, ...)
          })

setMethod("GeneSetCollection",
          signature=signature(
            object="ExpressionSet",
            idType="missing",
            setType="KEGGCollection"),
          function(object, ..., idType, setType) {
              idType <- AnnotationIdentifier(annotation(object))
              genes <- as.list(.getAnnMap(idType, "PATH2PROBE"))
              genes <- .GSC_filter_by_probe(genes, featureNames(object))
              .GSC_KEGG_helper(genes, idType, setType, ...)
          })

.GSC_GO_helper <- function(genes, idType, setType, ...) {
    ## filter on evidence codes
    evidenceCode = evidenceCode(setType)
    eviOk <- lapply(lapply(genes, names),
                    "%in%", evidenceCode)
    genes <- mapply("[", genes, eviOk)
    ugenes <- lapply(genes, unique)
    ugenes <- ugenes[sapply(ugenes, length) != 0]

    organism <- .getAnnMap(idType, "ORGANISM")
    gss <- mapply(function(ids, setName, collectionType, ...) {
        GeneSet(ids,
                setName=setName,
                collectionType=GOCollection(
                  goIds=setName,
                  evidenceCode=evidenceCode(collectionType)),
                ...)
    }, ugenes, names(ugenes), MoreArgs=list(
                                collectionType=setType,
                                geneIdType=idType,
                                organism=organism, ...))
    GeneSetCollection(gss)
}

setMethod("GeneSetCollection",
          signature=signature(
            object="missing",
            idType="AnnotationIdentifier",
            setType="GOCollection"),
          function(object, ..., idType, setType) {
              genes <- as.list(.getAnnMap(idType, "GO2PROBE"))
              .GSC_GO_helper(genes, idType=idType, setType=setType, ...)
          })

setMethod("GeneSetCollection",
          signature=signature(
            object="ExpressionSet",
            idType="missing",
            setType="GOCollection"),
          function(object, ..., idType, setType) {
              idType <- AnnotationIdentifier(annotation(object))
              genes <- as.list(.getAnnMap(idType, "GO2PROBE"))
              genes <- .GSC_filter_by_probe(genes, featureNames(object))
              .GSC_GO_helper(genes, idType=idType, setType=setType, ...)
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

## incidence

setMethod("incidence",
          signature=signature(
            x="GeneSetCollection"),
          function(x, ...) {
              args <- c(x, ...)
              .incidence(lapply(args, geneIds),
                         lapply(args, setName))
          })
