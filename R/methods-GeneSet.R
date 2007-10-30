.constructors_GeneSet("GeneSet", required=character(0))

## Rationale: 'initialize' is always called with .Object, constructed
## from the object prototype. The documentation for 'new' indicates
## that an un-named argument will be 'from [a] class[es] that this
## class extends', with the intention that unnamed arguments serve as
## templates to over-ride prototypes, but be over-written by specific
## named (slot) arguments.
## 
## Partial matching means that a named argument before ... without a
## corresponding named argument in a call to 'new' will pick up the
## template. To circumvent this partial matching, all named (slot)
## arguments in initialize should come _after_ ...
## 
## Template slots are meant to over-ride default slots provided in the
## prototype.  Named arguments requiring defaults need to capture
## these values from the template, rather than from the
## prototype. Providing a second argument before ..., called
## .Template, captures the (first) template argument, and named
## arguments after ... use this for defaults, rather than
## .Object. Since the template may not always be provided, .Template
## is in turn supplied with a default, i.e., .Object.
## 
## There are two additional outcomes of this scheme.
## 
## It is possible to use 'new' to update multiple slots, new('Obj',
## <obj>, slot1=<>, slot2=<>). This can be an efficient way to update
## many slots simultaneously, since 'new' is now relatively efficient
## in slot assignments (avoiding most unnecessary copying).
## 
## This also provides a way to perform 'atomic' updates, so that it is
## not necessary to construct invalid objects in the process of
## complex transformations.
## 
## An alternative approach would use named arguments to capture any
## incoming arguments before a .Object <- callNextMethod(.Object,
## ...), followed by slot assignment to .Object. This implies copying
## for each assignment, and precludes certain types of validity
## checking (e.g., when a named argument is meant to be 'required')

setMethod("initialize",
          signature=signature(.Object="GeneSet"),
          function(.Object, .Template=.Object, ...,
                   ## additional args, manipulated by method
                   setIdentifier=.uniqueIdentifier(),
                   setName=.Template@setName,
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

.GETTERS_GeneSet <- c("geneIdType", "geneIds", "setIdentifier",
                      "setName", description="shortDescription",
                      "longDescription", "organism", "pubMedIds", "urls",
                      "contributor", setVersion="version",
                      "creationDate", "collectionType")

.getters("GeneSet", .GETTERS_GeneSet)

.SETTERS_GeneSet <-
    .GETTERS_GeneSet["geneIdType" != .GETTERS_GeneSet]

.GeneSet_setters("GeneSet", .SETTERS_GeneSet[.SETTERS_GeneSet!="setIdentifier"])

.setters("GeneSet", "setIdentifier")    # no need to set unique identifier!

## convert between GeneIdentifier types

setReplaceMethod("geneIdType",
                 signature=signature(
                   object="GeneSet",
                   value="character"),
                 function(object, verbose=FALSE, value) {
                     tag <- tryCatch({
                         do.call(value, list())
                     }, error=function(err) {
                         stop(sprintf("could not create geneIdType tag of '%s'",
                                      value))
                     })
                     mapIdentifiers(object, tag, geneIdType(object), verbose=verbose)
                 })

setReplaceMethod("geneIdType",
                 signature=signature(
                   object="GeneSet",
                   value="GeneIdentifierType"),
                 function(object, verbose=FALSE, value) {
                     mapIdentifiers(object, value, geneIdType(object),
                                    verbose=verbose)
                 })

## subset

setMethod("[",
          signature=signature(
            x="GeneSet", i="numeric"),
          function(x, i, j, ..., drop=TRUE) {
              if (any(duplicated(i)))
                  stop("duplicate index: ",
                       paste(i[duplicated(i)], collapse=" "))
              geneIds <- geneIds(x)[i]
              if (any(is.na(geneIds)))
                  stop("unmatched index: ",
                       paste(i[is.na(geneIds)], collapse=" "))
              geneIds(x) <- geneIds
              x
          })

setMethod("[",
          signature=signature(
            x="GeneSet", i="character"),
          function(x, i, j, ..., drop=TRUE) {
              idx <- pmatch(i, geneIds(x))
              if (any(is.na(idx)))
                  stop(sprintf("unmatched / duplicate geneIds: '%s'",
                               paste(i[is.na(idx)], collapse="', '")))
              geneIds(x) <- geneIds(x)[idx]
              x
          })

setMethod("[[",
          signature=signature(
            x="GeneSet", i="numeric"),
          function(x, i, j, ...) {
              geneIds(x)[[i]]
          })

setMethod("[[",
          signature=signature(
            x="GeneSet", i="character"),
          function(x, i, j, ...) {
              idx <- pmatch(i, geneIds(x))
              if (is.na(idx))
                  stop("unmatched gene: ", i)
              geneIds(x)[[idx]]
          })

setMethod("$",
          signature=signature(x="GeneSet"),
          function(x, name) {
              i <- pmatch(name, geneIds(x), duplicates.ok=FALSE)
              if (is.na(i))
                  stop("unmatched gene: ", i)
              geneIds(x)[i]
          })

## Logic operations

.checkGeneSetLogicTypes <- function(x, y, functionName) {
    tx <- geneIdType(x)
    ty <- geneIdType(y)
    if (!(is(tx, class(ty)) || is(ty, class(tx))))
        stop(functionName, " incompatible GeneSet geneIdTypes;",
             "\n\tgot: ", class(tx), ", ", class(ty))
}

.geneSetIntersect <- function(x, y) {
   .checkGeneSetLogicTypes(x, y, "'&' or 'intersect'")
    new(class(x), x,
        setName=.glue(setName(x), setName(y), " & "),
        urls=.unique(urls(x), urls(y)),
        geneIds=intersect(geneIds(x), geneIds(y)),
        collectionType=intersect(collectionType(x), collectionType(y)))
}

.geneSetUnion <- function(x, y) {
    .checkGeneSetLogicTypes(x, y, "'|' or 'union'")
    new(class(x), x,
        setName=.glue(setName(x), setName(y), " | "),
        urls = .unique(urls(x), urls(y)),
        geneIds=union(geneIds(x), geneIds(y)),
        collectionType=union(collectionType(x), collectionType(y)))
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
              geneIds <- intersect(geneIds(e1), e2)
              new(class(e1), e1,
                  setName=.glue(setName(e1), "<character>", " & "),
                  geneIds=geneIds)
          })

setMethod("|",
          signature=signature(e1="GeneSet", e2="GeneSet"),
          function(e1, e2) .geneSetUnion(e1, e2))

setMethod("|",
          signature=signature(e1="GeneSet", e2="character"),
          function(e1, e2) {
              geneIds <- union(geneIds(e1), e2)
              new(class(e1), e1,
                  setName=.glue(setName(e1), "<character>", " | "),
                  geneIds=geneIds)
          })

setMethod("Logic",
          signature=signature(e1="character", e2="GeneSet"),
          function(e1, e2) callGeneric(e2, e1))

setMethod("setdiff",
          signature=signature(x="GeneSet", y="GeneSet"),
          function(x, y) {
              .checkGeneSetLogicTypes(x, y, "'setdiff'")
              geneIds=setdiff(geneIds(x), geneIds(y))
              new(class(x), x,
                  setName=.glue(setName(x), setName(y), " - "),
                  geneIds=setdiff(geneIds(x), geneIds(y)),
                  collectionType=setdiff(collectionType(x), collectionType(y)))
          })

## incidence

.incidence <- function(gidList, gnmList) {
    uids <- unique(unlist(gidList))
    isIn <- lapply(gidList,
                   function(g, u) match(u, g, nomatch=0),
                   uids)
    t(matrix(as.integer(unlist(isIn) > 0),
             ncol=length(gidList),
             dimnames=list(uids, unlist(gnmList))))
}

setMethod("incidence",
          signature=signature(
            x="GeneSet"),
          function(x, ...) {
              args <- c(x, ...)
              .incidence(lapply(args, geneIds),
                         lapply(args, setName))
          })

## show / description

.showGeneSet <- function(object) {
    cat("setName:", setName(object), "\n")
    cat("geneIds:",
        paste(selectSome(geneIds(object), maxToShow=4),
              collapse=", "),
        paste("(total: ", length(geneIds(object)), ")\n",
              sep=""),
        sep=" ")
    show(geneIdType(object))
    show(collectionType(object))
}

setMethod("show",
          signature=signature(object="GeneSet"),
          function(object) {
              .showGeneSet(object)
              cat("details: use 'details(object)'\n")
          })

setMethod("details",
          signature=signature(object="GeneSet"),
          function(object) {
              .showGeneSet(object)
              cat("setIdentifier: ", setIdentifier(object), "\n", sep="")
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
                  "setVersion: ", as(setVersion(object), "character"), "\n",
                  "creationDate: ", creationDate(object), "\n",
                  sep="")
          })
