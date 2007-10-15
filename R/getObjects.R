## file: source
## node: XPATH node identifier (syntax: http://www.w3.org/TR/xpath#path-abbrev)
## handler: function of 1 argument converting each node to GeneSet
.fromXML <- function(file, node, handler, ...) {
    res <- xmlTreeParse(file, useInternalNodes=TRUE, ...)
    geneSets <- getNodeSet(res, node, fun=handler)
    free(res)
    geneSets
}

.toXML <- function(geneSet, handler, con, ...) {
    saveXML(handler(geneSet), file=con, ...)
}

## Broad gene set, see http://www.broad.mit.edu/gsea/
.BROAD_SEPARATOR = ","                  # not specficied in DTD

.BroadXMLNodeToGeneSet_factory <- function(file) {
    ## state
    .mkSplit <- function(x) {
        if (is.null(x)) character(0)
        else unlist(strsplit(x, .BROAD_SEPARATOR))
    }
    url <- NULL
    if (length(file)==1) {
        isUri <- grep("^(http|ftp|file)://", file)
        if (length(isUri)==1 && isUri==1)
            url <- file
        else if (file.exists(file)) {
            full <- path.expand(file)
            if (length(grep("^/", full))==1)
                url <- paste("file:/", full, sep="")
            else if (.Platform$OS.type=="windows")
                url <- paste("file:///", full, sep="")
            else
                url = full
        }
    }
    ## handler: XMLNode -> GeneSet
    symbolId <- SymbolIdentifier()
    function(node) {
        attrs <- as.list(xmlAttrs(node))
        args <- list(symbolId,
                     setName=attrs[["STANDARD_NAME"]],
                     setIdentifier=attrs[["SYSTEMATIC_NAME"]],
                     geneIds=.mkSplit(attrs[["MEMBERS_SYMBOLIZED"]]),
                     organism=attrs[["ORGANISM"]],
                     urls= c(getBroadSet=url,
                       .mkSplit(attrs[["EXTERNAL_DETAILS_URL"]])),
                     collectionType={
                         categories <-
                             .mkSplit(attrs[["CATEGORY_CODE"]])
                         category <- subcategory <- as.character(NA)
                         if (length(categories)>=1)
                             category <- categories[[1]]
                         if (length(categories)>=2)
                             subcategory <- categories[[2]]
                         if (length(categories)>2)
                             warning("Broad 'CATEGORY_CODE' too long; using first two elements")
                         new("BroadCollection",
                             category=mkScalar(category),
                             subCategory=mkScalar(subcategory))
                     },
                     contributor=attrs[["CONTRIBUTOR"]],
                     pubMedIds=attrs[["PMID"]],
                     shortDescription=attrs[["DESCRIPTION_BRIEF"]],
                     longDescription=attrs[["DESCRIPTION_FULL"]],
                     TAGS=NULL, MESH=NULL, CHIP=NULL, MEMBERS=NULL)
        args <- args[!sapply(args, is.null)]
        do.call("GeneSet", args)
    }
}

.GeneSetToBroadXMLNode <- function(geneSet) {
    ## return text representation
    xmlNode("GENESET",
            attrs=c(
              STANDARD_NAME=setName(geneSet),
              SYSTEMATIC_NAME=setIdentifier(geneSet),
              ORGANISM=organism(geneSet),
              EXTERNAL_DETAILS_URL=paste(urls(geneSet),
                collapse=.BROAD_SEPARATOR),
              CATEGORY_CODE={
                  ct <- collectionType(geneSet)
                  category <- bcCategory(ct)
                  subcategory <- bcSubCategory(ct)
                  paste(if (is.na(category)) NULL else category,
                        if (is.na(subcategory)) NULL else subcategory,
                        collapse=.BROAD_SEPARATOR, sep="")
              },
              CONTRIBUTOR=contributor(geneSet),
              PMID=pubMedIds(geneSet),
              DESCRIPTION_FULL=longDescription(geneSet),
              DESCRIPTION_BRIEF=description(geneSet),
              MEMBERS_SYMBOLIZED=paste(geneIds(geneSet), collapse=",")
##               TAGS="",
##               MESH="",
##               CHIP="",
##               MEMBERS=""
              ))
}

asBroadUri <- function(name,
                       base="http://www.broad.mit.edu/gsea/msigdb/cards") {
    paste(base, "/", name, ".xml", sep="")
}

getBroadSets <- function(uri, ...) {
    factories <- sapply(uri, .BroadXMLNodeToGeneSet_factory)
    tryCatch({
        geneSets <- unlist(mapply(.fromXML, uri, "//GENESET", factories,
                                  SIMPLIFY=FALSE, USE.NAMES=FALSE))
    }, error=function(err) {
        stop("'getBroadSets' failed to create gene sets:\n  ",
             conditionMessage(err),
             call.=FALSE)
    })
    GeneSetCollection(geneSets)
}

toBroadXML <- function(geneSet, con = stdout(), ...) {
    if (!is(collectionType(geneSet), "BroadCollection"))
        .stopf("toBroadXML requires 'BroadCollection', got '%s'",
               class(collectionType(geneSet)))
    tryCatch({
        .toXML(geneSet, .GeneSetToBroadXMLNode, con=con, ...)
    }, error=function(err) {
        .stopf("toBroadXML failed to create XML: %s\n",
               conditionMessage(err))
    })
}

## 
##  OBO collections
## 

## Read 'obo' files for GO ids
.fromOBO <- function(src) {
    ## Parse OBO into 'stanza' and 'kv' (key-value) tables.
    ## VERY NAIVE
    parser <- list(blank_line="^\\s*$",
                   comment_line="^\\s*!",
                   comment="\\s*!.*$",
                   stanza="^\\[(.+)\\]",
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
    list(stanza=stanza, kv=kv)
}

.idsFromOBO <- function(stanza) {
    row.names(stanza)[stanza$Value=="Term"]
 }

getOBOCollection <- function(uri, ...) {
    res <- .fromOBO(uri)
    OBOCollection(.idsFromOBO(res$stanza), ...)
}
