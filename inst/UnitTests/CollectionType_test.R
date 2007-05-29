test_CollectionType_Constructors <- function() {
    ## do they exist and return the correct class?
    constructors <-
        GSEABase:::.nameAll(GSEABase:::.CONSTRUCTORS_CollectionType)
    for (i in seq_along(constructors)) {
        res <- do.call(constructors[[i]], list())
        checkTrue(validObject(res))
        checkTrue(is(res, constructors[[i]]))
    }
}
