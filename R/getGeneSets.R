## file: source
## node: XPATH node identifier (syntax: http://www.w3.org/TR/xpath#path-abbrev)
## handler: function of 1 argument converting each node to GeneSet
.fromXML <- function(file, node, handler, ...) {
    res <- xmlTreeParse(file, useInternalNodes=TRUE, ...)
    geneSets <- getNodeSet(res, node, fun=handler)
    free(res)
    GeneSetCollection(geneSets)
}

.toXML <- function(geneSet, handler, con, ...) {
    saveXML(handler(geneSet), file=con, ...)
}

## Broad gene set, see http://www.broad.mit.edu/gsea/
.BROAD_SEPARATOR = ","                  # not specficied in DTD

.BroadXMLNodeToGeneSet_factory <- function(file) {
    ## state
    .mkSplit <- function(x) unlist(strsplit(x, .BROAD_SEPARATOR))
    url <- NULL
    if (length(file)==1) {
        isUri <- grep("^(http|ftp|file)://", file)
        if (length(isUri)==1 && isUri==1)
            url <- file
        else if (file.exists(file)) {
            full <- path.expand(file)
            if (length(grep("^/", full))==1)
                url <- paste("file:/", full, sep="")
        }
    }
    ## handler: XMLNode -> GeneSet
    function(node) {
        attrs <- as.list(xmlAttrs(node))
        args <- list(new("SymbolIdentifier"),
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
                     TAGS=NULL,
                     MESH=NULL,
                     CHIP=NULL,
                     MEMBERS=NULL)
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

getBroadSets <- function(file, ...) {
    tryCatch({
        .fromXML(file, "//GENESET", .BroadXMLNodeToGeneSet_factory(file), ...)
    }, error=function(err) {
        stop("'getBroadSets' failed to create gene set:",
             "\n  ", conditionMessage(err),
             call.=FALSE)
    })
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

