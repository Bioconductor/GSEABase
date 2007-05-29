getClassOfSlot <- function(klass, slot) {
    as.character(getClass(klass)@slots[[slot]])
}

getterCheck <- function(getter, slotName, obj) {
    expectedClass <- getClassOfSlot(class(obj), slotName)
    checkTrue(is(getter(obj), expectedClass))
}

do_GeneSet_getter_check <- function(obj) {
    getters <- GSEABase:::.nameAll(GSEABase:::.GETTERS_GeneSet)
    for (g in names(getters)) {
        getterCheck(get(g), getters[[g]], obj)
    }
}

do_GeneSet_setter_check <- function(obj) {
    setters <- GSEABase:::.nameAll(GSEABase:::.SETTERS_GeneSet)
    for (s in names(setters)) {
        ss <- paste(s, "<-", sep="")
        obj <- do.call(ss,
                       list(obj, new(class(slot(obj, setters[[s]])))))
        checkTrue(validObject(obj,complete=TRUE))
    }
}

test_GS_MakeNoType <- function() {
    gs <- GeneSet(genes=letters[1:5],
                  setIdentifier="unique!",
                  setName="TestSet",
                  shortDescription="Test Gene Set No. 1",
                  longDescription="This is a gene set used for testing.",
                  organism="AlienX5.11",
                  pubMedIds=c("1", "2"),
                  urls=c("http://bioconductor.org"),
                  contributor="A.U. Thor")

    ## Basic accessor testing
    checkEquals(mkScalar("unique!"), setIdentifier(gs))
    checkEquals(mkScalar("TestSet"), setName(gs))

    do_GeneSet_getter_check(gs)
    do_GeneSet_setter_check(gs)
}

test_GS_MakeString <- function() {
    gs <- GeneSet("EntrezIdentifier",
                  genes=letters[1:5],
                  setIdentifier="unique!",
                  setName="TestSet",
                  shortDescription="Test Gene Set No. 1",
                  longDescription="This is a gene set used for testing.",
                  organism="AlienX5.11",
                  pubMedIds=c("1", "2"),
                  urls=c("http://bioconductor.org"),
                  contributor="A.U. Thor")
    checkEquals(mkScalar("unique!"), setIdentifier(gs))
    checkEquals(mkScalar("TestSet"), setName(gs))

    do_GeneSet_getter_check(gs)
    do_GeneSet_setter_check(gs)
}

test_GS_MakeType <- function() {
    gs <- GeneSet(new("EntrezIdentifier"),
                  genes=letters[1:5],
                  setIdentifier="unique!",
                  setName="TestSet",
                  shortDescription="Test Gene Set No. 1",
                  longDescription="This is a gene set used for testing.",
                  organism="AlienX5.11",
                  pubMedIds=c("1", "2"),
                  urls=c("http://bioconductor.org"),
                  contributor="A.U. Thor")

    checkEquals(mkScalar("unique!"), setIdentifier(gs))
    checkEquals(mkScalar("TestSet"), setName(gs))

    do_GeneSet_getter_check(gs)
    do_GeneSet_setter_check(gs)
}

test_GS_RequiredArgsToNew <- function() {
    checkException(GeneSet(new("EntrezIdentifier"),
                           genes=letters[1:5],
                           ## no setIdentifier
                           setName="TestSet",
                           shortDescription="Test Gene Set No. 1",
                           organism="AlienX5.11",
                           pubMedIds=c("1", "2"),
                           urls=c("http://bioconductor.org"),
                           contributor="A.U. Thor"),
                   silent=TRUE)
}

test_GS_MakeFromExpressionSet <- function() {
    data(sample.ExpressionSet)
    gs <- GeneSet(sample.ExpressionSet, setName="123",
                   setIdentifier="456")
    checkTrue(all(genes(gs)==featureNames(sample.ExpressionSet)))
    checkTrue(is(setType(gs), "AnnotationIdentifier"))
    checkTrue(is(collectionType(gs), "ExpressionSetCollection"))
    checkTrue(description(gs)==
              experimentData(sample.ExpressionSet)@title)
    checkTrue(longDescription(gs) ==
              abstract(experimentData(sample.ExpressionSet)))
    checkTrue(urls(gs) ==
              experimentData(sample.ExpressionSet)@url)
    checkTrue(contributor(gs) ==
              experimentData(sample.ExpressionSet)@name)
    do_GeneSet_getter_check(gs)
    do_GeneSet_setter_check(gs)
}

test_GS_setdiffExport <- function() {
    checkIdentical(environment(setdiff),
                   environment(GSEABase::setdiff))
}

test_GS_LogicalNonOverlapping <- function() {
    ## non-overlapping
    gss <- getBroadSets(system.file("extdata", "Broad.xml",
                                   package="GSEABase"))
    gs1 <- gss[[1]]
    gs2 <- gss[[2]]
    gs12 <- gs1 & gs2
    checkTrue(length(genes(gs12))==0)
    gs12 <- gs1 | gs2
    checkTrue(length(genes(gs12))==length(c(genes(gs1), genes(gs2))))
    checkTrue(all(genes(gs1) %in% genes(gs12)))
    checkTrue(all(genes(gs2) %in% genes(gs12)))
    checkIdentical(genes(GSEABase::setdiff(gs12, gs1)), genes(gs2))
    checkIdentical(genes(GSEABase::setdiff(gs12, gs2)), genes(gs1))
}

test_GS_LogicalOverlapping <- function() {
    gss <- getBroadSets(system.file("extdata", "Broad.xml",
                                   package="GSEABase"))
    gs1 <- gss[[1]]
    gs2 <- GeneSet(type=setType(gs1),
                   setName="123", setIdentifier="456",
                   genes=c(
                     sample(genes(gs1), 20),
                     sample(genes(gss[[2]]), 20)))
    checkTrue(all(genes(gs1 | gs2) %in%
                  genes(GSEABase::setdiff(gs1, gs2) |
                        (gs1 & gs2) |
                        GSEABase::setdiff(gs2, gs1))))
}
