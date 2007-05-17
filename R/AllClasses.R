## GeneIdentifierType

setClass("GeneIdentifierType", contains="VIRTUAL",
         representation=representation(
           type = "ScalarCharacter"))

setClass("Untyped",
         contains = "GeneIdentifierType",
         prototype = prototype(
           type = new("ScalarCharacter", "Untyped")))

setClass("EntrezIdentifier",
         contains = "GeneIdentifierType",
         prototype = prototype(
           type = new("ScalarCharacter", "Entrez")))

setClass("SymbolIdentifier",
         contains = "GeneIdentifierType",
         prototype = prototype(
           type=new("ScalarCharacter", "Symbol")))

setClass("PFAMIdentifier",
         contains = "GeneIdentifierType",
         prototype = prototype(
           type = new("ScalarCharacter", "PFAM")))

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
           uniqueIdentifier = "ScalarCharacter",
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
           type = new("Untyped"),
           version = new("Versions", "0.0.1"),
           collectionType = new("AdHocCollection")))
