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
        !(category %in% c("c1", "c2", "c3", "c4", "c5", "c6", "c7", "h")))
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
                         c4="c4 (Computational)",
                         c5="c5 (GO)",
                         c6="c6 (Oncogenic Pathway Activation Modules)",
                         c7="c7 (Immunologic Signatures)"),
                         h="h (Hallmark)", "\n",
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

.checkCode <- function(term, codes) {
    codeOk <- term %in% codes
    if (!all(codeOk))
        .stopf("%s invalid: %s\n  valid codes: %s",
               deparse(substitute(term)),
               paste(term[!codeOk], collapse=", "),
               paste(codes, collapse=", "))
    if ("ANY" %in% term)
      term <- codes[!codes %in% c("ANY", NA)]
    term
}

.checkGOEvidenceCodes <- function(evidenceCode) {
    codes <- c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "ISS", "ISO",
               "ISA", "ISM", "IGC", "IBA", "IBD", "IKR", "IRD", "RCA",
               "TAS", "NAS", "IC", "ND", "IEA", "ANY", NA)
    .checkCode(evidenceCode, codes)
}

.checkGOOntologyCode <- function(ontology) {
    .checkCode(ontology, c("CC", "MF", "BP", "ANY", NA))
}

.GOFilterIds <- function(ids, ontology, err=FALSE) {
    terms <- mget(ids, getAnnMap("TERM", "GO"), ifnotfound=NA_character_)
    termsOk <- !is.na(terms)
    if (err && !all(termsOk))
      .stopf("GO ids not found: '%s'",
             paste(ids[!termsOk], collapse="' '"))
    terms <- terms[termsOk]
    ids <- ids[termsOk]
    idsOk <- unlist(lapply(terms, Ontology)) %in% ontology
    if (err && !all(idsOk))
      .stopf("ontology '%s' does not contain ids '%s'",
             paste(ontology, collapse="' '"),
             paste(ids[!idsOk], collapse="' '"))
    ids[idsOk]
}

GOCollection <- function(ids=character(0),
                         evidenceCode="ANY",
                         ontology="ANY", ..., err=FALSE) {
    evidenceCode <- .checkGOEvidenceCodes(evidenceCode)
    ontology <- .checkGOOntologyCode(ontology)
    ids <- .GOFilterIds(ids, ontology, err=err)
    new("GOCollection", ids=ids, ...,
        evidenceCode=evidenceCode, ontology=ontology)
}

.SETTERS_GOCollection <- .GETTERS_GOCollection <-
    c("evidenceCode", "ontology")
.getters("GOCollection", .GETTERS_GOCollection)

setMethod("[",                          # e.g., x[evidenceCode="TAS"]
          signature=signature(
            x="GOCollection",
            i="missing",
            j="missing"),
          function(x, i, j, ...,
                   ids=GSEABase::ids(x),
                   evidenceCode=GSEABase::evidenceCode(x),
                   ontology=GSEABase::ontology(x),
                   drop=TRUE) {
              if (!missing(evidenceCode))
                  .collection_subset_chk("evidenceCode", evidenceCode,
                                         evidenceCode(x))
              if (!missing(ontology))
                  ids <- .GOFilterIds(ids, ontology, err=FALSE)
              callNextMethod(x,,,...,
                             ids=ids,
                             evidenceCode=evidenceCode,
                             ontology=ontology,
                             drop=drop)
          })

setMethod("show",
          signature=signature(object="GOCollection"),
          function(object) {
              callNextMethod()
              cat("  evidenceCode:", evidenceCode(object), "\n")
              cat("  ontology:", ontology(object), "\n")
          })
