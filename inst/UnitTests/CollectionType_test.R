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
    idx <- 200:220
    map <- getAnnMap("PFAM", annotation(eset))
    tbl <- toTable(map[featureNames(eset)[idx]])
    tbl <- tbl[!is.na(tbl$PfamId),]

    gsc <- GeneSetCollection(eset[idx,], setType=PfamCollection())
    checkIdentical(length(unique(tbl$PfamId)), length(gsc))
    x <- with(tbl, lapply(split(probe_id, PfamId), unique))
    checkIdentical(x[sort(names(x))],
                   geneIds(gsc)[order(names(gsc))])
}

test_CollectionType_PrositeCollection_ESet <-
    function()
{
    data(sample.ExpressionSet)
    idx <- 200:220
    eset <- sample.ExpressionSet[idx,]
    map <- getAnnMap("PROSITE", annotation(eset))
    tbl <- toTable(map[featureNames(eset)])
    tbl <- tbl[!is.na(tbl$ipi_id),]

    gsc <- GeneSetCollection(eset, setType=PrositeCollection())

    checkIdentical(length(unique(tbl$ipi_id)), length(gsc))
    x <- with(tbl, lapply(split(probe_id, ipi_id), unique))
    checkIdentical(x[sort(names(x))],
                   geneIds(gsc)[order(names(gsc))])
}

test_CollectionType_ChrlocCollection_ESet <-
    function()
{
    data(sample.ExpressionSet)
    eset <- sample.ExpressionSet
    idx <- 200:220
    map <- getAnnMap("CHRLOC", annotation(eset))
    tbl <- toTable(map[featureNames(eset)[idx]])

    gsc <- GeneSetCollection(eset[idx,], setType=ChrlocCollection())

    checkIdentical(nrow(unique(tbl[,2:3])), length(gsc))
    checkIdentical(length(unique(tbl[[1]])),
                   length(unique(unlist(geneIds(gsc)))))
    checkIdentical(sort(do.call(paste, c(unique(tbl[,3:2]), sep=":"))),
                   sort(names(gsc)))
}

test_CollectionType_PfamCollection_org <-
    function()
{
    idx <- 1:10
    map <- getAnnMap("PFAM", "org.Hs.eg")[idx]
    ids <- ls(map)
    tbl <- toTable(map)

    gsc <- GeneSetCollection(ids,
                             idType=AnnotationIdentifier("org.Hs.eg.db"),
                             setType=PfamCollection())
    len <- with(tbl, length(unique(PfamId[!is.na(PfamId)])))
    checkIdentical(len, length(gsc))
    gids <- unique(tbl[!is.na(tbl$PfamId),"gene_id"])
    x <- with(tbl, split(gene_id, PfamId))
    checkTrue(all(mapply(setequal, x, geneIds(gsc)[names(x)])))
}

test_CollectionType_PrositeCollection_org <-
    function()
{
    idx <- 1:10
    map <- getAnnMap("PROSITE", "org.Hs.eg")[idx]
    ids <- ls(map)
    tbl <- toTable(map)
    tbl <- tbl[!is.na(tbl[["ipi_id"]]),]

    gsc <- GeneSetCollection(ids,
                             idType=AnnotationIdentifier("org.Hs.eg.db"),
                             setType=PrositeCollection())

    checkIdentical(length(unique(tbl[["ipi_id"]])), length(gsc))
    x <- with(tbl, lapply(split(gene_id, ipi_id), unique))
    checkIdentical(x[order(names(x))],
                   geneIds(gsc)[order(names(gsc))])
}

test_CollectionType_ChrlocCollection_org <-
    function()
{
    idx <- 1:10
    map <- getAnnMap("CHRLOC", "org.Hs.eg")[idx]
    ids <- ls(map)
    tbl <- toTable(map)

    gsc <- GeneSetCollection(ids,
                             idType=AnnotationIdentifier("org.Hs.eg.db"),
                             setType=ChrlocCollection())
    checkIdentical(nrow(unique(tbl[,2:3])), length(gsc))
    grp <- do.call(paste, c(tbl[3:2], sep=":"))
    x <- with(tbl, lapply(split(gene_id, grp), unique))
    checkIdentical(x[sort(names(x))],
                   geneIds(gsc)[order(names(gsc))])
}
