.broadSets <- function()
    getBroadSets(system.file("extdata", "Broad.xml", package="GSEABase"))

.check_subset_Ok <- function(geneIds, gs) {
    checkTrue(validObject(gs))
    checkIdentical(length(geneIds), length(geneIds(gs)))
    checkTrue(all(geneIds==geneIds(gs)))
}

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
    gs <- GeneSet(geneIds=letters[1:5],
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
                  geneIds=letters[1:5],
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
                  geneIds=letters[1:5],
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
                           geneIds=letters[1:5],
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
    checkTrue("123" == setName(gs))
    checkTrue("456" == setIdentifier(gs))
    checkTrue(all(geneIds(gs)==featureNames(sample.ExpressionSet)))
    checkTrue(is(geneIdType(gs), "AnnotationIdentifier"))
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
    gss <- .broadSets()
    gs1 <- gss[[1]]
    gs2 <- gss[[2]]
    gs12 <- gs1 & gs2
    checkTrue(length(geneIds(gs12))==0)
    gs12 <- gs1 | gs2
    checkTrue(length(geneIds(gs12))==length(c(geneIds(gs1), geneIds(gs2))))
    checkTrue(all(geneIds(gs1) %in% geneIds(gs12)))
    checkTrue(all(geneIds(gs2) %in% geneIds(gs12)))
    checkIdentical(geneIds(GSEABase::setdiff(gs12, gs1)), geneIds(gs2))
    checkIdentical(geneIds(GSEABase::setdiff(gs12, gs2)), geneIds(gs1))
}

test_GS_LogicalOverlapping <- function() {
    gss <- getBroadSets(system.file("extdata", "Broad.xml",
                                   package="GSEABase"))
    gs1 <- gss[[1]]
    gs2 <- GeneSet(type=geneIdType(gs1),
                   setName="123", setIdentifier="456",
                   geneIds=c(
                     sample(geneIds(gs1), 20),
                     sample(geneIds(gss[[2]]), 20)))
    checkTrue(all(geneIds(gs1 | gs2) %in%
                  geneIds(GSEABase::setdiff(gs1, gs2) |
                        (gs1 & gs2) |
                        GSEABase::setdiff(gs2, gs1))))
}

test_GS_subset <- function() {
    gs <- .broadSets()[[1]]

    geneIds <- geneIds(gs)[4:1]
    .check_subset_Ok(geneIds, gs[4:1])
    .check_subset_Ok(geneIds, gs[ geneIds ])

    max <- length(geneIds(gs))
    checkException(gs[(max-1):(max+1)], silent=TRUE)
    checkException(gs["adfas"], silent=TRUE)
}

test_GS_subset2 <- function() {
    gs <- .broadSets()[[1]]

    geneIds <- geneIds(gs)[[2]]
    checkTrue(geneIds == gs[[2]])
    checkTrue(geneIds == gs[[ geneIds ]])

    max <- length(geneIds(gs))
    checkException(gs[["sdfsf"]], silent=TRUE)
    checkException(gs[[ max+1 ]], silent=TRUE)

    checkTrue(geneIds == do.call("$",list(gs, geneIds)))
}
