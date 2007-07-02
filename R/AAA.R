
.checkRequired <- function(required, provided) {
    idx <- which(!(required %in% provided))
    if (length(idx) > 0)
        stop("missing required argument(s): '",
             paste(required[idx], collapse="', '"), "'")
}

## for a vector c(CONSTRUCTOR=CLASS, ...) create a
## functionCONSTRUCTOR(...) calling new(CLASS, ...). Missing
## CONSTRUCTOR are filled with CLASS
.constructors_Simple <- function(klasses, required=NULL) {
    klassnames <- names(.nameAll(klasses))
    args <- .nameAll(c(required, "...")) # convenience of automatic matching
    iargs <- sapply(args, function(y) alist(y=)$y) # input args as pairlist
    oargs <- sapply(args, as.symbol)    # output args
    for (cl in seq_along(klasses))
        eval(substitute({
            f <- function() {
                .checkRequired(REQUIRED, names(match.call()))
                do.call("new", c(CLASS, OARGS))
            }
            formals(f) <- IARGS
            assign(CONSTRUCTOR, f, envir=topenv())
        }, list(CONSTRUCTOR = klassnames[[cl]],
                CLASS = klasses[[cl]],
                IARGS=iargs,
                OARGS=oargs,
                REQUIRED=required)))
}

## constructors for GeneSet and derived classes, with required fields.
.constructors_GeneSet<- function(klass, required) {
    ## construct the arg list of symbols with no defaults
    ## constructor input arguments: type, name, ...
    args <- .nameAll(c("type", required, "..."))
    iargs <- sapply(args, function(y) alist(y=)$y) # input args as pairlist
    oargs <- sapply(args[-1], as.symbol) # output args
    eval(substitute({
        if (!isGeneric(CLASS))
            setGeneric(CLASS,
                       function(type, ...) standardGeneric(CLASS))
        
        f <- function() {
            .checkRequired(REQUIRED, names(match.call()))
            do.call("new", c(CLASS, OARGS))
        }
        formals(f) <- IARGS
        setMethod(CLASS, signature = signature(type = "missing"), f)

        f <- function() {
            .checkRequired(REQUIRED, names(match.call()))
            do.call("new", c(CLASS, list(geneIds=type), OARGS))
        }
        formals(f) <- IARGS
        setMethod(CLASS, signature = signature(type = "character"), f)

        f <- function(){
            .checkRequired(REQUIRED, names(match.call()))
            do.call("new", c(CLASS, geneIdType=type, OARGS))
        }
        formals(f) <- IARGS
        setMethod(CLASS,
                  signature=signature(type="GeneIdentifierType"), f)

        f <- function(type, ...) {
            .checkRequired(REQUIRED, names(match.call()))
            organism <- 
                tryCatch({
                    pkg <- annotation(type)
                    if (length(pkg) == 1 && nchar(pkg) > 0 &&
                        .requireQ(pkg))
                        get(paste(pkg, "ORGANISM", sep=""))
                    else
                        ""
                }, error=function(err) "")
            new(CLASS,
                geneIdType = new("AnnotationIdentifier",
                  annotation = annotation(type)),
                geneIds = featureNames(type),
                shortDescription = experimentData(type)@title,
                longDescription = abstract(type),
                organism = organism,
                pubMedIds = pubMedIds(experimentData(type)),
                urls = experimentData(type)@url,
                contributor = experimentData(type)@name,
                collectionType = ExpressionSetCollection(),
                ...)
        }
        formals(f) <- IARGS
        setMethod(CLASS, signature = signature(type="ExpressionSet"), f)
    }, list(CLASS = klass, REQUIRED=required, IARGS=iargs, OARGS=oargs)))
}

.getters <- function(klass, slots) {
    slots <- .nameAll(slots)
    ## standard getters. 'where' default is topenv(parent.frame())
    ## which on package load is the package name space
    for (i in seq(along=slots)) {
        eval(substitute({
            if (!isGeneric(GENERIC))
                setGeneric(GENERIC,
                           function(object) standardGeneric(GENERIC))
            setMethod(GENERIC,
                      signature=signature(object=CLASS),
                      function(object) slot(object, SLOT))
        }, list(CLASS = klass,
                GENERIC = names(slots)[[i]],
                SLOT = slots[[i]])))
    }
}

.setters <- function(klass, slots) {
    slots <- .nameAll(slots)
    for (i in seq(along=slots)) {
        eval(substitute({

            if (!isGeneric(SETTER))
                setGeneric(SETTER, function(object, value)
                           standardGeneric(SETTER))
            if (getSlots(CLASS)[[SLOT]] == "ScalarCharacter")
                setReplaceMethod(GENERIC,
                                 signature=signature(
                                   object=CLASS,
                                   value="character"),
                                 function(object, value) {
                                     slot(object, SLOT) <- mkScalar(value)
                                     object
                                 })
            else
                setReplaceMethod(GENERIC,
                                 signature=signature(
                                   object=CLASS,
                                   value=getSlots(CLASS)[[SLOT]]),
                                 function(object, value) {
                                     slot(object, SLOT) <- value
                                     object
                                 })
        }, list(CLASS=klass,
                GENERIC=names(slots)[[i]],
                SETTER=paste(names(slots)[[i]], "<-", sep=""),
                SLOT=slots[[i]])))
    }
}
