## Constructors

OBOCollection <- function(ids=character(0),
                          evidenceCode="ANY",
                          ontology="ANY", ...) {
    evidenceCode <- .checkGOEvidenceCodes(evidenceCode)
    ontology <- .checkGOOntologyCode(ontology)
    if (any(duplicated(ids)))
        .stopf("OBO 'ids' contains duplicates: %s",
               paste(ids[duplicated(ids)], collapse=" "))
    .stanza <- data.frame(ids,
                          value=rep("Term", length(ids)),
                          row.names="ids")
    goMap <- getAnnMap("TERM", "GO")
    value <- rep(NA_character_, length(ids))
    goTerms <- mget(ids, goMap, ifnotfound=NA)
    ok <- !is.na(goTerms)
    if (sum(ok)>0)
        value[ok] <-  sapply(goTerms[ok], Term)
    .kv <- data.frame(stanza_id=ids,
                      key=rep("name", length(ids)),
                      value=value)
    new("OBOCollection", .stanza=.stanza, .kv=.kv,
        ids=ids, evidenceCode=evidenceCode, ontology=ontology, ...)
}

setValidity("OBOCollection",
            function(object) {
                msg <- NULL
                ids <- ids(object)
                stanza <- .stanza(object)
                okIds <- row.names(stanza)[stanza$value=="Term"]
                if (!all(ids %in% okIds)) {
                    bad <- ids[!ids %in% okIds]
                    str <- sprintf("identifiers in 'ids' but not '.stanza': %s",
                                   paste(bad, collapse=" "))
                    msg <- validMsg(msg, str)
                }
                if (!all(okIds %in% ids)) {
                    bad <- okIds[!okIds %in% ids]
                    str <- sprintf("identifiers in '.stanza' but not 'id': %s",
                                   paste(bad, collapse=" "))
                    msg <- validMsg(msg, str)
                }
                if (is.null(msg)) TRUE else msg
            })

## accessors, NOT exported
.SETTERS_OBOCollection <- .GETTERS_OBOCollection <-
    c(".stanza", ".kv", .obo_subset=".subset")
.setters("OBOCollection", .SETTERS_OBOCollection)
.getters("OBOCollection", .GETTERS_OBOCollection)

setMethod("subsets",
          signature=signature(
            object="OBOCollection"),
          function(object, display="named") {
              subsets <- .obo_subset(object)
              k <- row.names(subsets)
              v <- subsets[["value"]]
              switch(display,
                     "named"={
                         names(v) <- k
                         v
                     },
                     "key"=k,
                     "value"=v,
                     "full"={
                         if (length(k)==0)
                             character(0)
                         else
                             paste(k, " (", v, ")", sep="")
                     },
                     stop("'subsets' display must be missing, 'named', 'key', 'value', or 'full'"))
          })

setMethod("[",
          signature=signature(
            x="OBOCollection",
            i="character",
            j="missing"),
          function(x, i, j, ..., drop=TRUE) {
              i <- unique(i)
              if (!all(i %in% subsets(x, "key"))) {
                  bad <- i[!i %in% subsets(x, "key")]
                  .stopf("'[' indicies must be valid subsets.\n  bad values: %s\n possible values: %s\n",
                         bad, paste(subsets(x, "key"), collapse=", ", sep=""))
              }
              stanza <- .stanza(x)
              subset <- .obo_subset(x)
              kv <- .kv(x)

              kv_idx <- kv$key=="subset" & kv$value %in% i
              ids <- unique(kv$stanza_id[kv_idx])
              stanza_idx <-
                  row.names(stanza) %in% ids | stanza$value != "Term"

              stanza <- stanza[stanza_idx,, drop=FALSE]
              subset <- subset[row.names(subset) %in% i,,drop=FALSE]
              kv <- kv[kv$stanza_id %in% row.names(stanza),,drop=FALSE]

              callNextMethod(x,,, ...,
                             .stanza=stanza, .subset=subset, .kv=kv,
                             ids=.OBOids(stanza))
          })

## Read 'obo' files for GO ids
.fromOBO <- function(src, ...) {
    ## Parse OBO into 'stanza' and 'kv' (key-value) tables.
    ## VERY NAIVE
    parser <- list(blank_line="^\\s*$",
                   comment_line="^\\s*!",
                   comment="\\s*!.*$",
                   stanza="^\\[(.+)\\]",
                   subsetdef="^(\\w+)\\s*\"(.*)\"",
                   kv="^([^:]*):\\s*(.*)")
    data <- readLines(src)
    comments <- grep(parser$comment_line, data)
    if (length(comments) > 0)
        data <- data[-comments]
    data <- sub(parser$comment, "", data)
    ## stanzas: parser$blank_line followed by parser$stanza
    stnz <- grep(parser$stanza, data)
    stnz <- stnz[stnz %in% (grep(parser$blank_line, data) + 1)]
    stanza <- data.frame(id=c(0,stnz),
                         value=c("Root",
                           sub(parser$stanza, "\\1", data[stnz])),
                         stringsAsFactors=FALSE)
    ## key-value pairs, mapped to stanza
    kv_id <- grep(parser$kv, data)
    stanza_id <- sapply(kv_id, function(x) {
        idx <- x > stanza$id
        stanza$id[xor(idx, c(idx[-1], FALSE))]
    })

    kv_pairs <- data[kv_id]
    kv_key=sub(parser$kv, "\\1", kv_pairs)
    kv_value=sub(parser$kv, "\\2", kv_pairs)

    ## a lot of work to map stanza row.names as corresponding id's
    id_keys <- kv_key=="id"
    stanza_idx <- match(stanza_id[id_keys], stanza$id)
    row.names(stanza)[c(1, stanza_idx)] <-
        c(".__Root__", kv_value[id_keys])
    stanza_id <- row.names(stanza)[match(stanza_id, stanza$id)]
    stanza$id <- NULL

    ## kv data frame
    kv <- data.frame(id=kv_id[!id_keys],
                     stanza_id=stanza_id[!id_keys],
                     key=kv_key[!id_keys],
                     value=kv_value[!id_keys],
                     stringsAsFactors=FALSE,
                     row.names="id")

    ## subset
    idx <- kv$stanza_id==".__Root__" & kv$key=="subsetdef"
    if (sum(idx)>0) {
        v <- kv$value[idx]
        subset_id <- sub(parser$subsetdef, "\\1", v)
        subset_value <- sub(parser$subsetdef, "\\2", v)
        kv <- kv[!idx,]
    } else {
        subset_id <- subset_value <- character(0)
    }
    subset <- data.frame(id=subset_id,
                         value=subset_value,
                         row.names="id",
                         stringsAsFactors=FALSE)
    new("OBOCollection",
        .stanza=stanza, .subset=subset, .kv=kv,
        ids=.OBOids(stanza),
        ...)
}

.OBOids <- function(stanza)
    row.names(stanza)[stanza$value=="Term"]

getOBOCollection <- function(uri, evidenceCode="ANY", ...) {
    evidenceCode <- .checkGOEvidenceCodes(evidenceCode)
    .fromOBO(uri, evidenceCode=evidenceCode, ...)
}

setAs("OBOCollection", "graphNEL",
      function(from) {
          if (!.requireQ("graph")) {
              stop("package 'graph' required but not available")
          }
          s <- .stanza(from)
          k <- .kv(from)
          nodes <- .OBOids(s)
          edgeL <- rep(list(list(edges=integer())), length(nodes))
          names(edgeL) <- nodes
          df <- k[k$key=="is_a" & k$value %in% nodes,
                  !names(k)=="key",drop=FALSE]
          f <- function(x) list(edges=as.character(x))
          sid <- as.factor(df$stanza_id)
          edgeL[levels(sid)] <- lapply(split(df$value, sid), f)
          new("graphNEL", nodes, edgeL, "directed")
      })

setAs("graphNEL", "OBOCollection",
      function(from) {
          ids <- nodes(from)
          obo <- OBOCollection(ids)
          kv <- .kv(obo)

          iedge <- inEdges(from)
          nedge <- sapply(iedge, length)
          value <- mapply(rep, names(iedge), nedge, USE.NAMES=FALSE)
          if (length(value)){
              value <- unlist(value)
              stanza_id <- unlist(iedge, use.names=FALSE)
              kv <- rbind(.kv(obo),
                          data.frame(stanza_id, key="is_a", value))
          }
          new("OBOCollection", obo, .kv=kv)
      })

setMethod("show",
          signature=signature(
            object="OBOCollection"),
          function(object) {
              callNextMethod()
              subsets <- subsets(object, "full")
              str <- sprintf("subsets: %s (%d total)\n",
                             paste(selectSome(subsets, maxToShow=4),
                                   collapse=", "), length(subsets))
              writeLines(strwrap(str, indent=2, exdent=4))
          })

## Group GO ids into GO_ontology-specific GO_slim categories
.GO_slim <- function(ids, GO_slim, GO_ontology="MF", verbose=FALSE) {
    ## Get GO_slim terms, restricted to GO_ontology
    GOTERM <- getAnnMap("TERM", "GO")
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
                        MF=getAnnMap("MFOFFSPRING", "GO"),
                        BP=getAnnMap("BPOFFSPRING", "GO"),
                        CC=getAnnMap("CCOFFSPRING", "GO"),
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
          function(idSrc, slimCollection, ontology, evidenceCode="ANY", ...,
                   verbose=FALSE) {
              map <- .getAnnMap(idSrc, "GO")
              res <- mget(featureNames(sample.ExpressionSet), map)
              res <- res[!is.na(res)]
              gids <- unlist(lapply(res, subListExtract, "GOID"))
              ecode <- unlist(lapply(res, subListExtract, "Evidence"))
              if (evidenceCode != "ANY")
                  gids <- gids[ecode %in% evidenceCode]
              .GO_slim(unique(gids), ids(slimCollection), ontology)
          })
