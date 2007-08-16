.checkRequired <- function(required, provided) {
    idx <- which(!(required %in% provided))
    if (length(idx) > 0)
        stop("missing required argument(s): '",
             paste(required[idx], collapse="', '"), "'")
}

## for a vector c(CONSTRUCTOR=CLASS, ...) create a
## functionCONSTRUCTOR(...) calling new(CLASS, ...). Missing
## CONSTRUCTOR are filled with CLASS
.constructors_Simple <- function(klasses, required=NULL, where=topenv()) {
    klassnames <- names(.nameAll(klasses))
    args <- .nameAll(c(required, "...")) # convenience of automatic matching
    iargs <- sapply(args, function(y) alist(y=)$y) # input args as pairlist
    oargs <- sapply(args, as.symbol)    # output args
    for (cl in seq_along(klasses))
        eval(substitute({
            f <- function() {
                GSEABase:::.checkRequired(REQUIRED, names(match.call()))
                do.call("new", c(CLASS, OARGS))
            }
            formals(f) <- IARGS
            assign(CONSTRUCTOR, f, envir=WHERE)
        }, list(CONSTRUCTOR = klassnames[[cl]],
                CLASS = klasses[[cl]],
                IARGS=iargs,
                OARGS=oargs,
                REQUIRED=required,
                WHERE=where)))
}

## constructors for GeneSet and derived classes, with required fields.
.constructors_GeneSet<- function(klass, required) {
    ## construct the arg list of symbols with no defaults
    ## constructor input arguments: type, name, ...
    args <- .nameAll(c("type", required, "...", "setIdentifier"))
    iargs <- sapply(args, function(y) alist(y=)$y) # input args as pairlist
    iargs[["setIdentifier"]] <- quote(.uniqueIdentifier())
    oargs <- sapply(args[-1], as.symbol) # output args
    eval(substitute({
        if (!isGeneric(CLASS))
            setGeneric(CLASS,
                       function(type, ...,
                                setIdentifier=.uniqueIdentifier()) {
                           standardGeneric(CLASS)
                       },
                       signature=c("type"))

        ## missing
        f <- function() {
            .checkRequired(REQUIRED, names(match.call()))
            do.call("new", c(CLASS, OARGS))
        }
        formals(f) <- IARGS
        setMethod(CLASS, signature = signature(type = "missing"), f)

        ## character
        f <- function() {
            .checkRequired(REQUIRED, names(match.call()))
            do.call("new", c(CLASS, list(geneIds=type), OARGS))
        }
        formals(f) <- IARGS
        setMethod(CLASS, signature = signature(type = "character"), f)

        ## GeneIdentifierType
        f <- function(){
            .checkRequired(REQUIRED, names(match.call()))
            do.call("new", c(CLASS, geneIdType=type, OARGS))
        }
        formals(f) <- IARGS
        setMethod(CLASS,
                  signature=signature(type="GeneIdentifierType"), f)

        ## GOCollection
        f <- function() {
            .checkRequired(REQUIRED, names(match.call()))
            if (!.requireQ("GO"))
                stop("library(GO) required but not available")
            ## look for geneIdType; default is EntrezId
            idType <- geneIdType(geneIdType)
            lookup <- get(paste("GO", toupper(idType),sep=""))
            ids <- mget(goIds(type), lookup, ifnotfound=as.character(NA))
            ids <- lapply(ids,
                          function(x, codes) x[names(x) %in% codes],
                          evidenceCode(type))
            geneIds <- unique(unlist(ids, use.names=FALSE))
            do.call("new",
                    c(CLASS,
                      geneIdType=geneIdType,
                      collectionType=type,
                      list(geneIds=geneIds),
                      OARGS))
        }
        formals(f) <- c(IARGS[-length(IARGS)],
                        geneIdType=quote(EntrezIdentifier()),
                        IARGS[length(IARGS)])
        setMethod(CLASS,
                  signature=signature(type="GOCollection"), f)

        ## ExpressionSet
        f <- function() {
            .checkRequired(REQUIRED, names(match.call()))
            organism <- 
                tryCatch({
                    pkg <- annotation(type)
                    if (length(pkg) == 1 && nchar(pkg) > 0 &&
                        .requireQ(pkg))
                        .getAnnMap(type, "ORGANISM")
                    else
                        ""
                }, error=function(err) "")
            do.call("new",
                    c(CLASS,
                      geneIdType = AnnotationIdentifier(annotation(type)),
                      list(geneIds = featureNames(type)),
                      shortDescription = experimentData(type)@title,
                      longDescription = abstract(type),
                      organism = organism,
                      pubMedIds = pubMedIds(experimentData(type)),
                      urls = experimentData(type)@url,
                      contributor = experimentData(type)@name,
                      collectionType = ExpressionSetCollection(),
                      OARGS))
        }
        formals(f) <- IARGS
        setMethod(CLASS, signature = signature(type="ExpressionSet"), f)
    }, list(CLASS = klass, REQUIRED=required, IARGS=iargs, OARGS=oargs)))
}

.getters <- function(klass, slots, where=topenv(parent.frame()), ...) {
    slots <- .nameAll(slots)
    for (i in seq_along(slots)) {
        eval(substitute({
            if (!isGeneric(GENERIC, where=where))
                setGeneric(GENERIC,
                           function(object) standardGeneric(GENERIC),
                           where=WHERE)
            setMethod(GENERIC,
                      signature=signature(object=CLASS),
                      function(object) slot(object, SLOT),
                      where=WHERE)
        }, list(CLASS = klass,
                GENERIC = names(slots)[[i]],
                SLOT = slots[[i]],
                WHERE = where)))
    }
}

.setters <- function(klass, slots, where=topenv(parent.frame()), ...) {
    slots <- .nameAll(slots)
    for (i in seq(along=slots)) {
        eval(substitute({
            if (!isGeneric(SETTER, where=where))
                setGeneric(SETTER, function(object, value)
                           standardGeneric(SETTER),
                           where = WHERE)
            if (getSlots(CLASS)[[SLOT]] == "ScalarCharacter")
                setReplaceMethod(GENERIC,
                                 signature=signature(
                                   object=CLASS,
                                   value="character"),
                                 function(object, value) {
                                     slot(object, SLOT) <- mkScalar(value)
                                     validObject(object)
                                     object
                                 },
                                 where = WHERE)
            else
                setReplaceMethod(GENERIC,
                                 signature=signature(
                                   object=CLASS,
                                   value=getSlots(CLASS)[[SLOT]]),
                                 function(object, value) {
                                     slot(object, SLOT) <- value
                                     validObject(object)
                                     object
                                 },
                                 where = WHERE)
        }, list(CLASS=klass,
                GENERIC=names(slots)[[i]],
                SETTER=paste(names(slots)[[i]], "<-", sep=""),
                SLOT=slots[[i]],
                WHERE=where)))
    }
}
