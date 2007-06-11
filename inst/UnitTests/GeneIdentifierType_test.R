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
    checkTrue("hgu95av2"==annotation(AnnotationIdentifier("hgu95av2")))
}


test_GeneIdentifierType_setType <- function() {
    data(sample.ExpressionSet)

    gs <- GeneSet(sample.ExpressionSet[100:110],
                  setName="123", setIdentifier="456")
    suppressWarnings(setType(gs) <- "EntrezIdentifier")
    checkTrue(validObject(gs, complete=TRUE))
    checkTrue(is(setType(gs), "EntrezIdentifier"))

    gs <- GeneSet(sample.ExpressionSet[100:110],
                  setName="123", setIdentifier="456")
    suppressWarnings(setType(gs) <- EntrezIdentifier())
    checkTrue(validObject(gs, complete=TRUE))
    checkTrue(is(setType(gs), "EntrezIdentifier"))

    ## duplicate gene names exception
    gs <- GeneSet(sample.ExpressionSet[100:200],
                  setName="123", setIdentifier="456")
    opt <- options(warn=2)
    on.exit(options(opt))
    checkException(setType(gs) <- EntrezIdentifier(), silent=TRUE)
}

test_GeneIdentifierType_mapIdentifiers_toAnnotation <- function() {
    gss <- getBroadSets(system.file("extdata", "Broad.xml", package="GSEABase"))
    suppressWarnings({
        res <- mapIdentifiers(gss[[1]], AnnotationIdentifier("hgu95av2"))
    })
    detach("package:hgu95av2")
    checkTrue(validObject(res))
    checkEquals(41, length(genes(res)))
}

test_GeneIdentifierType_mapIdentifiers_nullAmbiguity <- function() {
    ## Original bug: 
    ##     1: Ambiguous method selection for "mapIdentifiers", target "GeneSet#AnnotationIdentifier#NullIdentifier" (the first of the signatures shown will be used)
    ##     GeneSet#AnnotationIdentifier#GeneIdentifierType
    ##     GeneSet#GeneIdentifierType#NullIdentifier
    opts <- options(warn=2)
    on.exit(options(opts))
    gs <- GeneSet(setName="123", setIdentifier="345")
    setType(gs) <- AnnotationIdentifier("xyz")
    checkTrue(validObject(gs))
}

test_GeneIdentifierType_mapIdentifiers_toAnnotation_via_Dbi <- function()  {
    gss <- getBroadSets(system.file("extdata", "Broad.xml", package="GSEABase"))
    suppressMessages(suppressWarnings({
        res <- mapIdentifiers(gss[[1]], AnnotationIdentifier("hgu95av2db"))
    }))
    detach("package:hgu95av2db")
    checkTrue(validObject(res))
    checkEquals(41, length(genes(res)))
}
