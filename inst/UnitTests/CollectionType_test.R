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

test_CollectionType_Logic <- function() {
    ops <- c("|", "&", intersect, union, setdiff)
    ## same CollectionType
    gs1 <- NullCollection()
    sapply(sapply(ops, do.call, list(gs1, gs1)),
           checkIdentical, NullCollection())
    ## different CollectionType
    gs2 <- ExpressionSetCollection()
    sapply(sapply(ops, do.call, list(gs1, gs2)),
           checkIdentical, ComputedCollection())
    ## CollectionIdType always 'Computed'
    gs3 <- KEGGCollection(ids=letters[1:3])
    sapply(sapply(ops, do.call, list(gs2, gs3)),
           checkIdentical, ComputedCollection())
    sapply(sapply(ops, do.call, list(gs3, gs3)),
           checkIdentical, ComputedCollection())
}
