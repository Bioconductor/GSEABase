.checkRequired <- function(required, ...) {
    idx <- which(!(required %in% names(match.call())[-1]))
    if (length(idx) > 0)
        stop("missing required argument(s) '",
             paste(required[idx], sep="', '"), "'")
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
                              if (length(pkg) == 1 &&
                                  nchar(pkg) > 0 &&
                                  require(pkg, quietly=TRUE, character.only=TRUE))
                                  get(paste(pkg, "ORGANISM", sep=""))
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

.nameAll <- function(x) {
    ## name = slot; missing name uses slot for getter name
    if (length(names(x)))
      names(x) <- ifelse(nchar(names(x)) == 0, x, names(x))
    else
      names(x) <- x
    x
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
