getClassOfSlot <- function(klass, slot) {
    as.character(getClass(klass)@slots[[slot]])
}

getterCheck <- function(getter, slotName, obj) {
    expectedClass <- getClassOfSlot(class(obj), slotName)
    checkTrue(is(getter(obj), expectedClass))
}

do_GeneSet_getter_check <- function(obj) {
    getters <- GSEABase:::.GETTERS_GeneSet
    names(getters) <- ifelse(nchar(names(getters)) == 0,
                             getters, names(getters))
    for (g in names(getters)) {
        getterCheck(get(g), getters[[g]], obj)
    }
}

test_MakeNoType <- function() {
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
}

test_MakeString <- function() {
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
}

test_MakeType <- function() {
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
}

test_RequiredArgsToNew <- function() {
    checkException(
                   GeneSet(new("EntrezIdentifier"),
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

test_MakeFromExpressionSet <- function() {
    data(sample.ExpressionSet)
    res <- GeneSet(sample.ExpressionSet, setName="123",
                   setIdentifier="456")
    checkTrue(all(genes(res)==featureNames(sample.ExpressionSet)))
    checkTrue(is(setType(res), "AnnotationIdentifier"))
    checkTrue(is(collectionType(res), "AdHocCollection"))
    checkTrue(description(res)==
              experimentData(sample.ExpressionSet)@title)
    checkTrue(longDescription(res) ==
              abstract(experimentData(sample.ExpressionSet)))
    checkTrue(urls(res) ==
              experimentData(sample.ExpressionSet)@url)
    checkTrue(contributor(res) ==
              experimentData(sample.ExpressionSet)@name)
}
