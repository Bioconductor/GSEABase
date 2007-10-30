## constructors / getters / setters in AllClasses.R

## CollectionType

.collectionTypeLogic <- function(x, y)
    if (class(x)!=class(y)) {
        ComputedCollection()
    } else x

setMethod("Logic",
          signature=signature(
            e1="CollectionType",
            e2="CollectionType"),
          function(e1, e2) .collectionTypeLogic(e1, e2))

setMethod("intersect",
          signature=signature(
            x="CollectionType",
            y="CollectionType"),
          .collectionTypeLogic)

setMethod("union",
          signature=signature(
            x="CollectionType",
            y="CollectionType"),
          .collectionTypeLogic)

setMethod("setdiff",
          signature=signature(
            x="CollectionType",
            y="CollectionType"),
          .collectionTypeLogic)

setMethod("show",
          signature=signature(object="CollectionType"),
          function(object) {
              cat("collectionType:", collectionType(object), "\n")
          })

## CollectionIdType

.collectionIdTypeLogic <- function(x, y)
    ComputedCollection()

setMethod("Logic",
         signature=signature(
           e1="CollectionIdType",
           e2="CollectionIdType"),
          function(e1, e2) .collectionIdTypeLogic(e1, e2))

setMethod("intersect",
          signature=signature(
            x="CollectionIdType",
            y="CollectionIdType"),
          .collectionIdTypeLogic)

setMethod("union",
          signature=signature(
            x="CollectionIdType",
            y="CollectionIdType"),
          .collectionIdTypeLogic)

setMethod("setdiff",
          signature=signature(
            x="CollectionIdType",
            y="CollectionIdType"),
          .collectionIdTypeLogic)

## BroadCollection

BroadCollection <- function(category="c1", subCategory=NA, ...) {
    if (length(category)!=1 ||
        !(category %in% c("c1", "c2", "c3", "c4")))
        stop(sprintf("invalid BroadCollection category: '%s'",
                     paste(category, collapse="', '")))
    new("BroadCollection",
        category=mkScalar(category),
        subCategory=mkScalar(as.character(subCategory)))
}

.GETTERS_BroadCollection <-
    c(bcCategory="category", bcSubCategory="subCategory")

.getters("BroadCollection", .GETTERS_BroadCollection)

setMethod("show",
          signature=signature(object="BroadCollection"),
          function(object) {
              cat("collectionType: ", collectionType(object), "\n",
                  "  bcCategory: ",
                  switch(bcCategory(object),
                         c1="c1 (Positional)",
                         c2="c2 (Curated)",
                         c3="c3 (Motif)",
                         c4="c4 (Computational)"), "\n",
                  "  bcSubCategory:  ", bcSubCategory(object), "\n", sep="")
          })

## CollectionIdType

.collection_subset_chk <- function(lbl, vals, ok) {
    if (!is.character(vals))
        .stopf("'%s' must be character, but is '%s'",
               lbl, class(vals),
               call.=FALSE)
    if (!all(vals %in% ok)) {
        bad <- vals[!vals %in% ok]
        .stopf("'%s' must subset current values\n  bad %s: %s",
               lbl, lbl, paste(bad, collapse=" "),
               call.=FALSE)
    }
}

setMethod("[",                          # e.g., x[evidenceCode="TAS"]
          signature=signature(
            x="CollectionIdType",
            i="missing",
            j="missing"),
          function(x, i, j, ..., ids=GSEABase::ids(x), drop=TRUE) {
              if (!missing(ids))
                  .collection_subset_chk("ids", ids, ids(x))
              new(class(x), x, ..., ids=ids)
          })

setMethod("show",
          signature = signature(
            object = "CollectionIdType"),
          function(object) {
              callNextMethod()
              ids <- ids(object)
              cat("  ids: ",
                  paste(selectSome(ids, maxToShow=4),
                        collapse=", "),
                  " (", length(ids), " total)\n", sep="")
          })

## GOCollection

.checkGOEvidenceCodes <- function(evidenceCode) {
    codes <- c("IMP", "IPI", "TAS", "ISS", "IDA", "NAS", "IEA", "IGI",
               "RCA", "IEP", "IC", "NR", "ND", "ANY", NA)
    codeOk <- evidenceCode %in% codes
    if (!all(codeOk))
        .stopf("evidenceCode invalid: %s\n  valid codes: %s",
               paste(evidenceCode[!codeOk], collapse=", "),
               paste(codes, collapse=", "))
    if ("ANY" %in% evidenceCode)
        evidenceCode <- codes[!codes %in% c("ANY", NA)]
    evidenceCode
}

GOCollection <- function(ids=character(0),
                         evidenceCode="ANY", ...) {
    evidenceCode <- .checkGOEvidenceCodes(evidenceCode)
    new("GOCollection", ids=ids, evidenceCode=evidenceCode)
}

.SETTERS_GOCollection <- .GETTERS_GOCollection <-
    c("evidenceCode")
.getters("GOCollection", .GETTERS_GOCollection)

setMethod("[",                          # e.g., x[evidenceCode="TAS"]
          signature=signature(
            x="GOCollection",
            i="missing",
            j="missing"),
          function(x, i, j, ...,
                   evidenceCode=GSEABase::evidenceCode(x),
                   drop=TRUE) {
              if (!missing(evidenceCode))
                  .collection_subset_chk("evidenceCode", evidenceCode,
                                         evidenceCode(x))
              callNextMethod(x,,,...,evidenceCode=evidenceCode,
                             drop=drop)
          })

setMethod("show",
          signature=signature(object="GOCollection"),
          function(object) {
              callNextMethod()
              cat("  evidenceCode:", evidenceCode(object), "\n")
          })
