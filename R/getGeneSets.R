## 'handlers' is a list of named functons. Function names correspond
## to XML elements. An additional function named getSets retrieves
## sets created after parsing the document
.getGeneSets <- function(file, handlers, ...) {
    res <- xmlTreeParse(file, handlers=handlers, ...)
    handlers$getSets()
}

## Broad gene set, see http://www.broad.mit.edu/gsea/
.broadXMLToGeneSet <- function(file) {
    ## separator not specified in DTD
    .mkSplit <- function(x) unlist(strsplit(x, ","))
    .mkNull <- function(x) NULL
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
    ## Format of list entries: XML_ATTRIBUTE = c("slot", converter)
    map <- 
        list(STANDARD_NAME=c("setName", quote(mkScalar)),
             SYSTEMATIC_NAME=c("setIdentifier", quote(mkScalar)),
             ORGANISM=c("organism", quote(mkScalar)),
             EXTERNAL_DETAILS_URL=c("urls",
               function(x) c(getBroadSets=url, .mkSplit(x))),
             CATEGORY_CODE=c("collectionType",
               function(x) {
                   categories <- .mkSplit(x)
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
               }),
             CONTRIBUTOR=c("contributor", quote(mkScalar)),
             PMID=c("pubMedIds", quote(.mkSplit)),
             GEOID=c(NA, quote(.mkNull)),
             ## 'DESCRIPTION' in DTD
             DESCRIPTION_FULL=c("longDescription", quote(mkScalar)),
             DESCRIPTION_BRIEF=c("shortDescription", quote(mkScalar)),
             TAGS=c(NA, quote(.mkNull)),
             MESH=c(NA, quote(.mkNull)),
             MEMBERS_SYMBOLIZED=c("genes", quote(.mkSplit)),
             ## FIXME: 'original' source ?
             CHIP=c(NA, quote(.mkNull)),
             MEMBERS=c(NA, quote(.mkNull)))
    sets <- list()
    list(getSets=function() sets,
         GENESET=function(x) {
             ## Parse GENESET to object of class GeneSet
             attrs <- xmlAttrs(x)
             vattrs <-
                 lapply(names(attrs),
                        function(x) eval(map[[x]][[2]])(attrs[[x]]))
             names(vattrs) <- sapply(names(attrs),
                                     function(x) map[[x]][[1]])
             vattrs <- vattrs[!is.na(names(vattrs))]
             sets <<- append(sets,
                             do.call("GeneSet",
                                     c(new("EntrezIdentifier"), vattrs)))
             NULL
         })
}

getBroadSets <- function(file, ...) {
    tryCatch({
        .getGeneSets(file, .broadXMLToGeneSet(file), ...)
    }, error=function(err) {
        stop("'broadGeneSet' failed to create gene set:",
             "\n  ", conditionMessage(err),
             call.=FALSE)
    })
}
