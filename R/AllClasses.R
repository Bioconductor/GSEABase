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
           type = new("ScalarCharacter", "Entrez")))

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

setClass("BroadCollection",
         contains = "CollectionType",
         representation = representation(
           category = "ScalarCharacter",
           subCategory = "ScalarCharacter"),
         prototype = prototype(
           type = new("ScalarCharacter", "Broad"),
           category = new("ScalarCharacter", "none"),
           subCategory = new("ScalarCharacter", "none")))

## GeneSet

setClass("GeneSet",
         representation = representation(
           ## Gene set representation
           type = "GeneIdentifierType",
           genes = "character",
           ## Descriptive metadata
           setIdentifier = "ScalarCharacter",
           setName = "ScalarCharacter",
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
           type = new("NullIdentifier"),
           version = new("Versions", "0.0.1"),
           collectionType = new("AdHocCollection")),
         validity = function(object) {
             if (any(duplicated(genes(object))))
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
                 any(clen != length(genes(object))))
                 msg <- c(msg,"gene and color lengths differ")
             if (!("factor" %in% class(geneColor(object))) ||
                 !("factor" %in% class(phenotypeColor(object))))
                 msg <- c(msg, "gene- and phenotypeColor must be 'factor'")
             if (!is.null(msg))
                 msg
             else TRUE
         })
