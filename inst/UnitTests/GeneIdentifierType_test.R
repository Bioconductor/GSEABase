test_GeneIdentifierType_Constructors <- function() {
    ## do they exist and return the correct class?
    constructors <-
        GSEABase:::.nameAll(GSEABase:::.CONSTRUCTORS_GeneIdentifierType)
    constructors <- constructors[constructors!="AnnotationIdentifier"]
    for (i in seq_along(constructors)) {
        res <- do.call(constructors[[i]], list())
        checkTrue(validObject(res))
        checkTrue(is(res, constructors[[i]]))
    }

    ## Required slot for AnnotationIdentifier
    res <- AnnotationIdentifier(annotation="123")
    checkTrue(validObject(res))
    checkTrue(is(res, "AnnotationIdentifier"))
    checkException(AnnotationIdentifier(), silent=TRUE) # required arg missing
    checkTrue("hgu95av2"==geneType(AnnotationIdentifier("hgu95av2")))

    ## Optional slot for NullIdentifier
    checkTrue("<Ad hoc>"==geneType(NullIdentifier()))
    checkTrue("myType"==geneType(NullIdentifier("myType")))
}


test_GeneIdentifierType_mapIdentifiers <- function() {
    data(sample.ExpressionSet)

    gs <- GeneSet(sample.ExpressionSet[100:110],
                  setName="123", setIdentifier="456")
    setType(gs) <- "EntrezIdentifier"
    checkTrue(validObject(gs, complete=TRUE))
    checkTrue(is(setType(gs), "EntrezIdentifier"))

    gs <- GeneSet(sample.ExpressionSet[100:110],
                  setName="123", setIdentifier="456")
    setType(gs) <- EntrezIdentifier()
    checkTrue(validObject(gs, complete=TRUE))
    checkTrue(is(setType(gs), "EntrezIdentifier"))

    ## duplicate gene names exception
    gs <- GeneSet(sample.ExpressionSet[100:200],
                  setName="123", setIdentifier="456")
    checkException(setType(gs) <- EntrezIdentifier(), silent=TRUE)
}
