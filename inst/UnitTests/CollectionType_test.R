library(org.Hs.eg.db)
idConstructors <- names(slot(getClass("CollectionIdType"), "subclasses"))
simpleConstructors <- local({
    idTypes <- names(slot(getClass("CollectionType"), "subclasses"))
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

test_CollectionType_PfamCollection_ESet <-
    function()
{
    data(sample.ExpressionSet)
    eset <- sample.ExpressionSet
    gsc <- GeneSetCollection(eset[200:220,], setType=PfamCollection())
    checkIdentical(18L, length(gsc))
    checkIdentical(15L, length(unique(unlist(geneIds(gsc)))))
}

test_CollectionType_PrositeCollection_ESet <-
    function()
{
    data(sample.ExpressionSet)
    eset <- sample.ExpressionSet
    gsc <- GeneSetCollection(eset[200:220,], setType=PrositeCollection())
    checkIdentical(97L, length(gsc))
    checkIdentical(16L, length(unique(unlist(geneIds(gsc)))))
    checkIdentical(c("IPI00002903", "IPI00003325", "IPI00008359",
                     "IPI00008524", "IPI00009984", "IPI00018230"),
                     head(unique(names(gsc))))
}

test_CollectionType_ChrlocCollection_ESet <-
    function()
{
    data(sample.ExpressionSet)
    eset <- sample.ExpressionSet
    gsc <- GeneSetCollection(eset[200:220,], setType=ChrlocCollection())
    checkIdentical(17L, length(gsc))
    checkIdentical(14L, length(unique(unlist(geneIds(gsc)))))
    checkIdentical(c("11:-46655207", "1:-153201397", "12:-51448205",
                     "12:55018929", "1:-25561326", "1:47674275",
                     "18:-38577189", "19:-47917633", "19:59777068",
                     "19:59796924", "21:-35081967", "21:-35115443",
                     "4:71261081", "5:133478300", "5:133479214",
                     "5:133479325", "8:-101784319" ), names(gsc))
}

test_CollectionType_PfamCollection_org <-
    function()
{
    ids <- ls(org.Hs.egPFAM[1:10])
    gsc <- GeneSetCollection(ids,
                             idType=AnnotationIdentifier("org.Hs.eg.db"),
                             setType=PfamCollection())
    checkIdentical(10L, length(gsc))
    checkIdentical(6L, length(unique(unlist(geneIds(gsc)))))
    checkIdentical(c("1000", "1", "10000", "10000", "10000", "10",
                     "100", "1000", "100008586", "1000"),
                     unlist(lapply(gsc, geneIds)))
}

test_CollectionType_PrositeCollection_org <-
    function()
{
    ids <- ls(org.Hs.egPFAM[1:10])
    gsc <- GeneSetCollection(ids,
                             idType=AnnotationIdentifier("org.Hs.eg.db"),
                             setType=PrositeCollection())
    checkIdentical(19L, length(gsc))
    checkIdentical(6L, length(unique(unlist(geneIds(gsc)))))
    checkIdentical(c("10", "100008586", "1", "10000", "10000", "1000",
                     "100", "100008586", "1", "1", "1000", "1",
                     "1000", "10", "1000", "1000", "100", "10000",
                     "10000"), unlist(lapply(gsc, geneIds)))
}

test_CollectionType_ChrlocCollection_org <-
    function()
{
    ids <- ls(org.Hs.egPFAM[1:100])
    gsc <- GeneSetCollection(ids,
                             idType=AnnotationIdentifier("org.Hs.eg.db"),
                             setType=ChrlocCollection())
    checkIdentical(12L, length(gsc))
    checkIdentical(9L, length(unique(unlist(geneIds(gsc)))))
    checkIdentical(c("10003", "10000", "10000", "10001", "10002",
                     "1000", "1", "100", "10", "100008586",
                     "100008586", "100008586"), unlist(lapply(gsc,
                     geneIds)))
}
