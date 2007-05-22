.constructors <- function(klass)
    eval(substitute({
        setMethod(CLASS,
                  signature = signature(type = "missing"),
                  function(type, ...) {
                      new(CLASS, type=new("NullIdentifier"), ...)
                  })
        setMethod(CLASS,
                  signature = signature(type = "character"),
                  function(type, ...) {
                      new(CLASS, type=new(type), ...)
                  })
        setMethod(CLASS,
                  signature = signature(type="GeneIdentifierType"),
                  function(type, ...) {
                      new(CLASS, type = type, ...)
                  })
        setMethod(CLASS,
                  signature = signature(type="ExpressionSet"),
                  function(type, ...) {
                      organism <- 
                          tryCatch({
                              pkg <- annotation(type)
                              if (require(pkg, quietly=TRUE, character.only=TRUE))
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
    }, list(CLASS = klass)))

.getters <- function(klass, slots) {
    ## getter name = slot; missing name uses slot for getter name
    if (length(names(slots)))
      names(slots) <- ifelse(nchar(names(slots)) == 0,
                             slots, names(slots))
    else
      names(slots) <- slots
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
