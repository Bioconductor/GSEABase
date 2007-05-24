.checkRequired <- function(required, ...) {
    idx <- which(!(required %in% names(match.call())[-1]))
    if (length(idx) > 0)
        stop("missing required argument(s) '",
             paste(required[idx], collapse="', '"), "'")
}

.constructors <- function(klass, required)
    eval(substitute({
        setMethod(CLASS,
                  signature = signature(type = "missing"),
                  function(type, ...) {
                      .checkRequired(REQUIRED, ...)
                      new(CLASS, type=new("NullIdentifier"), ...)
                  })
        setMethod(CLASS,
                  signature = signature(type = "character"),
                  function(type, ...) {
                      .checkRequired(REQUIRED, ...)
                      new(CLASS, type=new("NullIdentifier"), ...)
                  })
        setMethod(CLASS,
                  signature = signature(type="GeneIdentifierType"),
                  function(type, ...) {
                      .checkRequired(REQUIRED, ...)
                      new(CLASS, type = type, ...)
                  })
        setMethod(CLASS,
                  signature = signature(type="ExpressionSet"),
                  function(type, ...) {
                      .checkRequired(REQUIRED, ...)
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
                  })
    }, list(CLASS = klass, REQUIRED=required)))

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
