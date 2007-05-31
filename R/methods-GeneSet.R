.constructors_GeneSet("GeneSet",
                      required=c("setName", "setIdentifier"))

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

.SETTERS_GeneSet <-
    .GETTERS_GeneSet["setType" != names(.GETTERS_GeneSet)]

.setters("GeneSet", .SETTERS_GeneSet)

## convert between GeneIdentifier types

setReplaceMethod("setType",
                 signature=signature(
                   object="GeneSet",
                   value="character"),
                 function(object, value) {
                     tag <- tryCatch({
                         do.call(value, list())
                     }, error=function(err) {
                         stop(sprintf("could not create setType tag of '%s'",
                                      value))
                     })
                     mapIdentifiers(setType(object), tag, object)
                 })

setReplaceMethod("setType",
                 signature=signature(
                   object="GeneSet",
                   value="GeneIdentifierType"),
                 function(object, value) {
                     mapIdentifiers(setType(object), value, object)
                 })

## subset

setMethod("[",
          signature=signature(
            x="GeneSet", i="numeric"),
          function(x, i, j, ..., drop=TRUE) {
              if (any(duplicated(i)))
                  stop("duplicate index: ",
                       paste(i[duplicated(i)], collapse=" "))
              genes <- genes(x)[i]
              if (any(is.na(genes)))
                  stop("unmatched index: ",
                       paste(i[is.na(genes)], collapse=" "))
              genes(x) <- genes
              x
          })

setMethod("[",
          signature=signature(
            x="GeneSet", i="character"),
          function(x, i, j, ..., drop=TRUE) {
              idx <- pmatch(i, genes(x))
              if (any(is.na(idx)))
                  stop(sprintf("unmatched / duplicate genes: '%s'",
                               paste(i[is.na(idx)], collapse="', '")))
              genes(x) <- genes(x)[idx]
              x
          })

setMethod("[[",
          signature=signature(
            x="GeneSet", i="numeric"),
          function(x, i, j, ...) {
              genes(x)[[i]]
          })

setMethod("[[",
          signature=signature(
            x="GeneSet", i="character"),
          function(x, i, j, ...) {
              idx <- pmatch(i, genes(x))
              if (is.na(idx))
                  stop("unmatched gene: ", i)
              genes(x)[[idx]]
          })

setMethod("$",
          signature=signature(x="GeneSet"),
          function(x, name) {
              i <- pmatch(name, genes(x), duplicates.ok=FALSE)
              if (is.na(i))
                  stop("unmatched gene: ", i)
              genes(x)[i]
          })

## Logic operations

.checkGeneSetLogicTypes <- function(x, y, functionName) {
    tx <- setType(x)
    ty <- setType(y)
    if (!(is(tx, class(ty)) || is(ty, class(tx))))
        stop(functionName, " incompatible GeneSet setTypes;",
             "\n\tgot: ", class(tx), ", ", class(ty))
}

.geneSetIntersect <- function(x, y) {
   .checkGeneSetLogicTypes(x, y, "'&' or 'intersect'")
    new(class(x), x,
        setIdentifier=setIdentifier(x),
        setName=.glue(setName(x), setName(y), " & "),
        urls=.unique(urls(x), urls(y)),
        genes=intersect(genes(x), genes(y)))
}

.geneSetUnion <- function(x, y) {
    .checkGeneSetLogicTypes(x, y, "'|' or 'union'")
    new(class(x), x,
        setIdentifier=setIdentifier(x),
        setName=.glue(setName(x), setName(y), " | "),
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
                  setName=.glue(setName(e1), "<character>", " & "),
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
                  setName=.glue(setName(e1), "<character>", " | "),
                  genes=genes)
          })

setMethod("Logic",
          signature=signature(e1="character", e2="GeneSet"),
          function(e1, e2) callGeneric(e2, e1))

setMethod("setdiff",
          signature=signature(x="GeneSet", y="GeneSet"),
          function(x, y) {
              .checkGeneSetLogicTypes(x, y, "'setdiff'")
              genes=setdiff(genes(x), genes(y))
              new(class(x), x,
                  setIdentifier=setIdentifier(x),
                  setName=.glue(setName(x), setName(y), " - "),
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
