## constructors / getters / setters in AllClasses.R

setMethod("show",
          signature=signature(object="CollectionType"),
          function(object) {
              cat("collectionType:", collectionType(object), "\n")
          })

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

setMethod("show",
          signature=signature(object="GOCollection"),
          function(object) {
              callNextMethod()
              cat("  evidenceCode:", evidenceCode(object), "\n")
          })

## OBOCollection

OBOCollection <- function(ids=character(0),
                          evidenceCode="ANY", ...) {
    evidenceCode <- .checkGOEvidenceCodes(evidenceCode)
    if (any(duplicated(ids)))
        .stopf("OBO 'ids' contains duplicates: %s",
               paste(ids[duplicated(ids)], collapse=" "))
    .stanza <- data.frame(ids,
                          value=rep("Term", length(ids)),
                          row.names="ids")
    goMap <- getAnnMap("TERM", "GO")
    .kv <- data.frame(stanza_id=ids,
                      key=rep("name", length(ids)),
                      value=sapply(mget(ids, goMap), Term))
    new("OBOCollection", .stanza=.stanza, .kv=.kv,
        ids=ids, evidenceCode=evidenceCode)
}

## Group GO ids into GO_ontology-specific GO_slim categories
.GO_slim <- function(ids, GO_slim, GO_ontology="MF", verbose=FALSE) {
    require("AnnotationDbi")
    require("GO.db")

    ## Get GO_slim terms, restricted to GO_ontology
    terms <- mget(GO_slim, GOTERM, ifnotfound=NA)
    if (any(is.na(terms))) {
        if (verbose)
            warning("GO_slim ids not found: ",
                    paste(names(terms)[is.na(terms)], collapse=" "))
        terms <- terms[!is.na(terms)]
    }
    terms <- terms[sapply(terms, Ontology)==GO_ontology]
    GO_slim <- names(terms)

    ## Use GO_ontology to find the required offspring
    OFFSPRING <- switch(GO_ontology,
                        MF=GOMFOFFSPRING,
                        BP=GOBPOFFSPRING,
                        CC=GOCCOFFSPRING,
                        stop("GO_ontology must be 'MF', 'BP', or 'CC'"))
    ## Get the offspring of GO_slim
    slim <- mget(GO_slim, OFFSPRING, ifnotfound=NA)
    slim <- slim[!is.na(slim)]
    ## Reverse the relationship: 'offspring' become keys, 'parents'
    ## become values. Select the sampled offspring
    ids <- unique(ids)
    samp <- revmap(slim)[ids]
    samp <- samp[!sapply(samp, is.null)]
    ## Count occurences of each slim
    cnt <- table(unlist(samp))
    ## Adjust for sample ids matching slim ids
    idx <- table(ids[which(ids %in% names(slim))])
    idx_n<- names(idx)
    cnt[idx_n] <- idx + ifelse(is.na(cnt[idx_n]), 0, cnt[idx_n])

    ## Prepare a data frame for results
    df <- data.frame(Slim=names(terms),
                     Count=0L, Percent=0,
                     Term=sprintf("%.35s%s",
                       sapply(terms, Term, USE.NAMES=FALSE),
                       ifelse(nchar(sapply(terms, Term))>35, "...", "")),
                     row.names=1)
    ## add our counts
    df[names(cnt),c("Count", "Percent")] <- c(cnt, 100*cnt/sum(cnt))
    df[order(row.names(df)),]
}

setMethod("goSlim",
          signature=signature(
            idSrc="GOCollection",
            slimCollection="GOCollection"),
          function(idSrc, slimCollection, ontology, ..., verbose=FALSE) {
              .GO_slim(ids(idSrc), ids(slimCollection), ontology, verbose)
          })

setMethod("goSlim",
          signature=signature(
            idSrc="ExpressionSet",
            slimCollection="GOCollection"),
          function(idSrc, slimCollection, ontology, evidenceCode="ANY", ..., verbose=FALSE) {
              map <- .getAnnMap(idSrc, "GO")
              res <- mget(featureNames(sample.ExpressionSet), map)
              res <- res[!is.na(res)]
              gids <- unlist(lapply(res, subListExtract, "GOID"))
              ecode <- unlist(lapply(res, subListExtract, "Evidence"))
              if (evidenceCode != "ANY")
                  gids <- gids[ecode %in% evidenceCode]
              .GO_slim(unique(gids), ids(slimCollection), ontology)
          })
