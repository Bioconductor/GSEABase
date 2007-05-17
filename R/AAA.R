.constructors <- function(klass)
    eval(substitute({
        setGeneric(CLASS,
                   function(type, ..., creationDate = date())
                   standardGeneric(CLASS))
        setMethod(CLASS,
                  signature = signature(type = "missing"),
                  function(type, ..., creationDate) {
                      new(CLASS, type = new("NullIdentifier"), ...,
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
    if (length(names(slots)))
      names(slots) <- ifelse(nchar(names(slots)) == 0,
                             slots, names(slots))
    else
      names(slots) <- slots
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
