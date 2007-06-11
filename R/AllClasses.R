## GeneIdentifierType

setClass("GeneIdentifierType", contains="VIRTUAL",
         representation=representation(
           type = "ScalarCharacter"))

setClass("NullIdentifier",
         contains = "GeneIdentifierType",
         prototype = prototype(
           type = new("ScalarCharacter", "NullIdentifier")))

setClass("AnnotationIdentifier",
         contains = "GeneIdentifierType",
         representation = representation(
           annotation = "ScalarCharacter"), # Bioconductor annotation package
         prototype = prototype(
           type = new("ScalarCharacter", "Annotation")))

setClass("EntrezIdentifier",
         contains = "GeneIdentifierType",
         prototype = prototype(
           type = new("ScalarCharacter", "EntrezId")))

setClass("SymbolIdentifier",
         contains = "GeneIdentifierType",
         prototype = prototype(
           type=new("ScalarCharacter", "Symbol")))

setClass("PfamIdentifier",
         contains = "GeneIdentifierType",
         prototype = prototype(
           type = new("ScalarCharacter", "Pfam")))

## CollectionType

setClass("CollectionType",
         contains = "VIRTUAL",
         representation = representation(
           type = "ScalarCharacter"))

setClass("AdHocCollection",
         contains = "CollectionType",
         prototype = prototype(
           type = new("ScalarCharacter", "Ad hoc")))

setClass("ExpressionSetCollection",
         contains = "CollectionType",
         prototype = prototype(
           type = new("ScalarCharacter", "ExpressionSet")))

setClass("KEGGCollection",
         contains = "CollectionType",
         prototype = prototype(
           type = new("ScalarCharacter", "KEGG")))

setClass("GOCollection",
         contains = "CollectionType",
         prototype = prototype(
           type = new("ScalarCharacter", "GO")))

setClass("BroadCollection",
         contains = "CollectionType",
         representation = representation(
           category = "ScalarCharacter",
           subCategory = "ScalarCharacter"),
         prototype = prototype(
           type = new("ScalarCharacter", "Broad"),
           category = new("ScalarCharacter", "c1"),
           subCategory = new("ScalarCharacter", as.character(NA))))

## GeneSet

setClass("GeneSet",
         representation = representation(
           ## Gene set representation
           geneIdType = "GeneIdentifierType",
           geneIds = "character",
           ## Descriptive metadata
           setName = "ScalarCharacter",
           setIdentifier = "ScalarCharacter",
           shortDescription = "ScalarCharacter",
           longDescription = "ScalarCharacter",
           organism = "ScalarCharacter",
           pubMedIds = "character",
           urls = "character",
           contributor = "character",
           version = "Versions",
           creationDate = "character",
           collectionType = "CollectionType"),
         prototype = prototype(
           setName = new("ScalarCharacter", "<undefined>"),
           setIdentifier = new("ScalarCharacter", "<undefined>"),
           geneIdType = new("NullIdentifier"),
           version = new("Versions", "0.0.1"),
           collectionType = new("AdHocCollection")),
         validity = function(object) {
             if (any(duplicated(geneIds(object))))
                 "gene symbols must be unique"
             else
                 TRUE
         })

setClass("GeneColorSet",
         contains = "GeneSet",
         representation = representation(
           phenotype = "ScalarCharacter",
           geneColor = "factor",
           phenotypeColor = "factor"),
         validity = function(object) {
             msg <- NULL
             clen <- c(length(geneColor(object)),
                       length(phenotypeColor(object)))
             if (any(clen > 0) &
                 any(clen != length(geneIds(object))))
                 msg <- c(msg,"gene and color lengths differ")
             if (!("factor" %in% class(geneColor(object))) ||
                 !("factor" %in% class(phenotypeColor(object))))
                 msg <- c(msg, "gene- and phenotypeColor must be 'factor'")
             if (!is.null(msg))
                 msg
             else TRUE
         })

## GeneSetCollection

setClass("GeneSetCollection",
         contains="list",
         validity = function(object) {
             msg <- NULL
             if (!all(sapply(object, is, "GeneSet")))
                 msg <- c(msg, "members must all be 'GeneSet' classes")
             tryCatch({
                 if (any(duplicated(names(object))))
                     msg <- c(msg, "each setName must be distinct")
                 }, error=function(err) {
                     msg <<- c(msg, conditionMessage(err))
                 })
             if (!is.null(msg))
                 msg
             else
                 TRUE
         })
