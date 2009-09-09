
## Information; hgu95av2ORGANISM, hgu95av2CHRLENGTHS        
## 
## Primary identifier: hgu95av2ACCNUM
## 
## Gene identifiers: hgu95av2ENTREZID hgu95av2ENZYME hgu95av2GENENAME
## hgu95av2PFAM hgu95av2PROSITE hgu95av2REFSEQ hgu95av2SYMBOL
## hgu95av2UNIGENE
## 
## Collections: hgu95av2OMIM hgu95av2PMID hgu95av2GO hgu95av2CHR
## hgu95av2CHRLOC hgu95av2MAP hgu95av2PATH
## 
## Reverse maps: hgu95av2ENZYME2PROBE hgu95av2GO2ALLPROBES
## hgu95av2GO2PROBE hgu95av2PATH2PROBE hgu95av2PMID2PROBE
## 
## Quality control: hgu95av2QC hgu95av2MAPCOUNTS
## Deprecated: hgu95av2LOCUSID

## GeneIdentifierType
.IdentifierClasses <- function(where) {
    setSimpleClass <- function(class, contains) {
        classId <- paste(class, "Identifier", sep="")
        setClass(classId,
                 contains = contains,
                 prototype = prototype(
                   type = new("ScalarCharacter", class)),
                 where = where)
    }

    setClass("GeneIdentifierType",
             contains="VIRTUAL",
             representation=representation(
               type = "ScalarCharacter",
               annotation = "ScalarCharacter"),
             where=where)
    .getters("GeneIdentifierType", c(geneIdType="type"), where=where)
    .setters("GeneIdentifierType", "annotation", where=where)

    ## Straight derivation from GeneIdentifierType
    annoIdentifiers <- c("Null", "Enzyme", "Genename", "Refseq", "Symbol",
                         "Unigene", "ENSEMBL")
    for (class in annoIdentifiers)
      setSimpleClass(class, "GeneIdentifierType")
    ## More complicated derviation from GeneIdentifierType:
    setClass("AnnotationIdentifier",    # 'annotation' slot
             contains = c("GeneIdentifierType"),
             prototype = prototype(
               type = new("ScalarCharacter", "Annotation")),
             where = where)
    setClass("EntrezIdentifier",        # special 'type'
             contains = "GeneIdentifierType",
             prototype = prototype(
               type = new("ScalarCharacter", "EntrezId")),
             where = where)
    idTypes <- names(slot(getClass("GeneIdentifierType"), "subclasses"))
    .constructors_Simple(idTypes, optional="annotation", where=where)
 }

.IdentifierClasses(topenv())

## Special class of Identifier for GOAllFrames
setClass("GOAllFrameIdentifier",
         contains="GeneIdentifierType",
         representation=representation(organism="ScalarCharacter"))  


## CollectionType
.CollectionClasses <- function(where) {
    setSimpleCollection <- function(class, contains) {
        classCollection <- paste(class, "Collection", sep="")
        setClass(classCollection,
                 contains = contains,
                 prototype = prototype(
                   type = new("ScalarCharacter", class)),
                 where = where)
    }

    ## simple collections
    setClass("CollectionType",
             contains = "VIRTUAL",
             representation = representation(
               type = "ScalarCharacter"),
             where = where)

    simpleCollections <- c("Null", "ExpressionSet", "Computed")
    for (class in simpleCollections)
        setSimpleCollection(class, "CollectionType")

    ## collections with ids
    setClass("CollectionIdType",
             contains=c("CollectionType", "VIRTUAL"),
             representation = representation(
               ids = "character"),
             where = where)
    
    idCollections <- c("KEGG", "OMIM", "PMID", "Chr", "Chrloc",
                       "Map", "Pfam", "Prosite")
    for (class in idCollections)
        setSimpleCollection(class, "CollectionIdType")

    ## other collections
    setClass("GOCollection",
             contains = "CollectionIdType",
             representation=representation(
               evidenceCode="character",
               ontology="character"),
             prototype = prototype(
               type = new("ScalarCharacter", "GO"),
               evidenceCode = NA_character_,
               ontology = NA_character_),
             where = where)

    setClass("OBOCollection",
             contains="GOCollection",
             representation=representation(
               .stanza="data.frame",
               .subset="data.frame",
               .kv="data.frame"),
             prototype=prototype(
               type=new("ScalarCharacter", "OBO"),
               .stanza=data.frame(id=character(0),
                 value=character(0), row.names="id"),
               .subset=data.frame(id=character(0),
                 value=character(0), row.names="id"),
               .kv=data.frame(stanza_id=character(0),
                 key=character(0), value=character(0))
               ),
             where=where)

    setClass("BroadCollection",
             contains = "CollectionType",
             representation = representation(
               category = "ScalarCharacter",
               subCategory = "ScalarCharacter"),
             prototype = prototype(
               type = new("ScalarCharacter", "Broad"),
               category = new("ScalarCharacter", "c1"),
               subCategory = new("ScalarCharacter", as.character(NA))),
             where = where)

    ## constructors / getters / setters
    ## (GOCollection and BroadCollection in methods-CollectionType.R)
    ## (OBOCollection in methods-OBOCollection.R)
    simpleCollections <- paste(simpleCollections, "Collection", sep="")
    .constructors_Simple(simpleCollections, where=where)
    .getters("CollectionType", c(collectionType="type"), where=where)

    idCollections <- paste(idCollections, "Collection", sep="")
    .constructors_Simple(idCollections, optional="ids", where=where)
    .getters("CollectionIdType", c("ids"), where=where)
    .setters("CollectionIdType", c("ids"), where=where)
}

.CollectionClasses(topenv())

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
           setName = new("ScalarCharacter", NA),
           setIdentifier = new("ScalarCharacter", NA),
           geneIdType = new("NullIdentifier"),
           version = new("Versions", "0.0.1"),
           collectionType = new("NullCollection")),
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
