idConstructors <- names(getSubclasses(getClass("CollectionIdType")))
simpleConstructors <- local({
    idTypes <- names(getSubclasses(getClass("CollectionType")))
    idTypes[!idTypes %in% c("CollectionIdType", idConstructors)]
})
constructors <- c(simpleConstructors, idConstructors)

test_CollectionType_Constructors <- function() {
    ## do they exist and return the correct class?
    for (i in seq_along(constructors)) {
        res <- do.call(constructors[[i]], list())
        checkTrue(validObject(res))
        checkTrue(is(res, constructors[[i]]))
    }

    ## BroadCollection
    checkTrue("c2" == bcCategory(BroadCollection(category="c2")))
    checkTrue("yyy" == bcSubCategory(BroadCollection(subCategory="yyy")))
}

test_getOBOCollection <- function() {
    fl <- system.file("extdata", "goslim_plant.obo", package="GSEABase")
    obj <- getOBOCollection(fl, evidenceCode="TAS")
    checkEquals(106, length(ids(obj)))
    checkEquals("TAS", evidenceCode(obj))
}
