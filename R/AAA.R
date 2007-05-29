
.checkRequired <- function(required, provided) {
    idx <- which(!(required %in% provided))
    if (length(idx) > 0)
        stop("missing required argument(s): '",
             paste(required[idx], collapse="', '"), "'")
}

## for a vector c(CONSTRUCTOR=CLASS, ...) create a
## functionCONSTRUCTOR(...) calling new(CLASS, ...). Missing
## CONSTRUCTOR are filled with CLASS
.constructors_Simple <- function(klasses) {
    klassnames <- names(.nameAll(klasses))
    for (cl in seq_along(klasses))
        eval(substitute({
            assign(CONSTRUCTOR,
                   function(...) new(CLASS, ...),
                   envir=topenv(parent.frame()))
        }, list(CONSTRUCTOR = klassnames[[cl]],
                CLASS = klasses[[cl]])))
}

## constructors for GeneSet and derived classes, with required fields.
.constructors_GeneSet<- function(klass, required) {
    ## construct the arg list of symbols with no defaults
    ## constructor input arguments: type, name, ...
    iargs <- sapply(.nameAll(c("type", required, "...")),
                    function(y) alist(y=)$y)
    ## arguments as seen by 'new': name=name, ...
    oargs <- sapply(c(.nameAll(required), "..."), as.symbol)
    eval(substitute({
        if (!isGeneric(CLASS))
            setGeneric(CLASS,
                       function(type, ...) standardGeneric(CLASS))
        
        f <- function() {
            .checkRequired(REQUIRED, names(match.call()))
            do.call("new", c(CLASS, type=new("NullIdentifier"), OARGS))
        }
        formals(f) <- IARGS
        setMethod(CLASS, signature = signature(type = "missing"), f)

        f <- function() {
            .checkRequired(REQUIRED, names(match.call()))
            do.call("new", c(CLASS, type=new(type), OARGS))
        }
        formals(f) <- IARGS
        setMethod(CLASS, signature = signature(type = "character"), f)

        f <- function(){
            .checkRequired(REQUIRED, names(match.call()))
            do.call("new", c(CLASS, type=type, OARGS))
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
                type = new("AnnotationIdentifier",
                  annotation = annotation(type)),
                genes = featureNames(type),
                shortDescription = experimentData(type)@title,
                longDescription = abstract(type),
                organism = organism,
                pubMedIds = pubMedIds(experimentData(type)),
                urls = experimentData(type)@url,
                contributor = experimentData(type)@name,
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
