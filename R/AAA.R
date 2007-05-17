.constructors <- function(klass)
    eval(substitute({
        setGeneric(CLASS,
                   function(type, ..., creationDate = date())
                   standardGeneric(CLASS))
        setMethod(CLASS,
                  signature = signature(type = "missing"),
                  function(type, ..., creationDate) {
                      new(CLASS, type = new("Untyped"), ...,
                          creationDate = creationDate)
                  })
        setMethod(CLASS,
                  signature = signature(type = "character"),
                  function(type, ..., creationDate) {
                      new(CLASS, type = new(type), ...,
                          creationDate = creationDate)
                  })
        setMethod(CLASS,
                  signature = signature(type="GeneIdentifierType"),
                  function(type, ..., creationDate) {
                      new(CLASS, type = type, ...,
                          creationDate = creationDate)
                  })
    }, list(CLASS = klass)))

.getters <- function(klass, slots) {
    ## getter name = slot; missing name uses slot for getter name
    autoName <-
        if (length(names(slots))) nchar(names(slots))==0
        else rep(TRUE, length(slots))
    names(slots)[autoName] <- slots[autoName]
    ## standard getters. 'where' default is topenv(parent.frame())
    ## which on package load is the package name space
    for (i in seq(along=slots)) {
        eval(substitute({
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
