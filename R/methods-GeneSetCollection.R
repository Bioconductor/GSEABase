setMethod("GeneSetCollection",
          signature=signature(
            object="GeneSet"),
          function(object, ...) new("GeneSetCollection", list(object, ...)))

setMethod("GeneSetCollection",
          signature=signature(
            object="list"),
          function(object, ...) new("GeneSetCollection", object))

## setMethod("GeneSetCollection",
##           signature=signature(
##             object="ExpressionSet"),
##           function(object,
##                    geneIdType=NullIdentifier(),
##                    collectionType=NullCollection(), ...) {
##               if (is(collectionType, "NullCollection")) {
##                   GeneSetCollection(GeneSet(object, ...))
##               } else if (is(collectionType, "KEGGCollection")) {
##                   probe2path <- mget(featureNames(object),
##                                      revmap(hgu95av2PATH2PROBE),
##                                      ifnotfound=list(NULL))
##                   probe2path <- probe2path[!sapply(probe2path, is.null)]

##                   hasPath <- sapply(featureNames(object), exists, revmap(hgu95av2PATH2PROBE))
##                   hasPath <- hasPath[hasPath]

##                   browser()
##               } else {
##                   stop("not yet implemented")
##               }
##           })

setMethod("GeneSetCollection",
          signature=signature(
            object="KEGGCollection"),
          function(object, annotation, ...) {
              if (missing(annotation))
                  stop("'annotation' package required")
              require(annotation, character.only=TRUE)
              annEnv <- paste(annotation, "PATH2PROBE", sep="")
              gss <- eapply(get(annEnv), GeneSet,
                            collectionType=KEGGCollection(),
                            geneIdType=AnnotationIdentifier(annotation))
              gss <- mapply("setName<-", gss, names(gss), USE.NAMES=FALSE)
              GeneSetCollection(gss)
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
