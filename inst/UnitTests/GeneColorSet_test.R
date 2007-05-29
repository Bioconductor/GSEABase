.colorBroadSets <- function() {
    gss <- getBroadSets(system.file("extdata", "Broad.xml",
                                    package="GSEABase"))
    lapply(gss, function(gs) {
        gcs <- GeneColorSet(gs, phenotype="<undefined>")
        df <- coloring(gcs)
        df[,"geneColor"] <- as.factor(rep(c("P", "M"), length=nrow(df)))
        df[,"phenotypeColor"] <- as.factor(rep(LETTERS[1:3], length=nrow(df)))
        coloring(gcs) <- df
        gcs
    })
}

test_GCS_ConstructorNoColorSetArgs <- function() {
    checkException(GeneColorSet(setIdentifier="123",
                                setName="Set name"),
                   silent=TRUE)
}

test_GCS_ConstructorDefaultArgs <- function() {
    gs <- GeneColorSet(setIdentifier="123",
                       setName="Set name",
                       phenotype="A phenotype")
    checkTrue(validObject(gs, complete=TRUE))
    checkTrue(length(genes(gs))==0)
    checkTrue(length(geneColor(gs))==0)
    checkTrue(length(phenotypeColor(gs))==0)
}

test_GCS_ConstructorAllColorArgs <- function() {
    ## appropriate default colors
    gs <- GeneColorSet(setIdentifier="123",
                       setName="Set name",
                       phenotype="A phenotype",
                       genes=LETTERS[1:24])
    checkTrue(validObject(gs, complete=TRUE))
    checkIdentical(genes(gs), LETTERS[1:24])
    checkTrue(length(geneColor(gs))==24)
    checkTrue(length(phenotypeColor(gs))==24)

    ## correct color contents
    gfactor <- factor(rep(c("high", "low"), 12))
    pfactor <- factor(rep(c("big", "medium", "small"), 8))
    gs <- GeneColorSet(setIdentifier="123",
                       setName="Set name",
                       phenotype="A phenotype",
                       genes=LETTERS[1:24],
                       geneColor=gfactor,
                       phenotypeColor=pfactor)
    checkTrue(validObject(gs, complete=TRUE))
    checkIdentical(geneColor(gs), gfactor)
    checkIdentical(phenotypeColor(gs), pfactor)
}

test_GCS_show <- function() {
    gs <- GeneColorSet(setIdentifier="123",
                       setName="Set name",
                       phenotype="A phenotype",
                       genes=LETTERS[1:24])
    con <- textConnection("tmp", open="w", local=TRUE)
    sink(con)
    on.exit(sink())
    show(gs)
}

test_GCS_colorizeReplace <- function() {
    gcs <- .colorBroadSets()[[1]]
    df <- coloring(gcs)
    gc <- as.factor(rep(c("Up", "Down"), length=nrow(df)))
    pc <- as.factor(rep(c("Bigger", "Smaller", "Same"),
                        length=nrow(df)))
    df$geneColor <- gc
    df$phenotypeColor <- pc
    coloring(gcs) <- df

    checkIdentical(geneColor(gcs), gc)
    checkIdentical(phenotypeColor(gcs), pc)
}

test_GCS_colorizeReplaceRetainGeneOrder <- function() {
    gcs <- .colorBroadSets()[[1]]
    ogenes <- genes(gcs)
    coloring(gcs) <- coloring(gcs)[sample(ogenes, length(ogenes)),]
    checkIdentical(genes(gcs), ogenes)
}

test_GCS_intersect <- function() {
    gcss <- .colorBroadSets()

    res <- GSEABase::intersect(gcss[[1]], gcss[[2]])
    checkTrue(validObject(res, complete=TRUE))
    checkTrue(length(genes(res))==0)
    checkTrue(length(urls(res))==4)
    checkIdentical(levels(geneColor(gcss[[1]])),
                   levels(geneColor(res)))
    checkIdentical(levels(phenotypeColor(gcss[[1]])),
                   levels(phenotypeColor(res)))

    res <- GSEABase::intersect(gcss[[1]], gcss[[1]])
    checkTrue(validObject(res, complete=TRUE))
    checkIdentical(genes(gcss[[1]]), genes(res))
    checkIdentical(urls(gcss[[1]]), urls(res))
    checkIdentical(levels(geneColor(gcss[[1]])),
                   levels(geneColor(res)))
    checkIdentical(levels(phenotypeColor(gcss[[1]])),
                   levels(phenotypeColor(res)))
}

test_GCS_intersectDifferentColors <- function() {
    gcs1 <- .colorBroadSets()[[1]]
    gcs2 <- gcs1
    phenotype(gcs2) <- paste(phenotype(gcs2), "A")
    geneColor(gcs2) <- factor(rep(c("Q", "R"), length=length(genes(gcs2))))
    ## warning about synthetic phenotype
    oldOpts <- options(warn=2)
    on.exit(options(oldOpts))
    checkException(res <- GSEABase::intersect(gcs1, gcs2), silent=TRUE)
    options(oldOpts)
    ## 
    suppressWarnings(res <- GSEABase::intersect(gcs1, gcs2))
    checkTrue(phenotype(res) != phenotype(gcs1))
    checkIdentical(2L, length(levels(geneColor(res))))
    checkIdentical(3L, length(levels(phenotypeColor(res))))
    checkTrue(!any(levels(geneColor(res)) == levels(geneColor(gcs2))))
}

test_GCS_union <- function() {
    gcss <- .colorBroadSets()

    res <- GSEABase::union(gcss[[1]], gcss[[2]])
    checkTrue(validObject(res, complete=TRUE))
    checkIdentical(length(genes(res)),
                   sum(sapply(gcss,
                              function(x) length(genes(x)))))
    checkTrue(all(genes(res) %in% c(genes(gcss[[1]]),
                                    genes(gcss[[2]]))))
    checkIdentical(levels(geneColor(gcss[[1]])),
                   levels(geneColor(res)))
    checkIdentical(levels(phenotypeColor(gcss[[1]])),
                   levels(phenotypeColor(res)))
    checkTrue(all(urls(res) %in% unlist(sapply(gcss, urls))))

    res <- GSEABase::union(gcss[[1]], gcss[[1]])
    checkTrue(validObject(res))
    checkIdentical(genes(res), genes(gcss[[1]]))
    checkIdentical(geneColor(res), geneColor(gcss[[1]]))
    checkIdentical(phenotypeColor(res), phenotypeColor(gcss[[1]]))
    checkIdentical(urls(gcss[[1]]), urls(res))
}

test_GCS_setdiff <- function() {
    gcss <- .colorBroadSets()

    res <- GSEABase::setdiff(gcss[[1]], gcss[[2]])
    checkTrue(validObject(res, complete=TRUE))
    checkIdentical(genes(gcss[[1]]), genes(res))
    checkIdentical(geneColor(gcss[[1]]), geneColor(res))
    checkIdentical(phenotypeColor(gcss[[1]]), phenotypeColor(res))
    
    res <- GSEABase::setdiff(gcss[[2]], gcss[[1]])
    checkTrue(validObject(res, complete=TRUE))
    checkIdentical(genes(gcss[[2]]), genes(res))
    checkIdentical(geneColor(gcss[[2]]), geneColor(res))
    checkIdentical(phenotypeColor(gcss[[2]]), phenotypeColor(res))

    res <- GSEABase::setdiff(gcss[[1]], gcss[[1]])
    checkTrue(validObject(res, complete=TRUE))
    checkTrue(length(genes(res))==0)
    checkIdentical(levels(geneColor(gcss[[1]])), levels(geneColor(res)))
    checkIdentical(levels(phenotypeColor(gcss[[1]])), levels(phenotypeColor(res)))
}

test_GCS_LogicalNonOverlapping <- function() {
    gcss <- .colorBroadSets()
    gs1 <- gcss[[1]]
    gs2 <- gcss[[2]]

    gs12 <- gs1 & gs2
    checkTrue(length(genes(gs12))==0)

    gs12 <- gs1 | gs2
    checkTrue(length(genes(gs12))==length(c(genes(gs1), genes(gs2))))
    checkTrue(all(genes(gs1) %in% genes(gs12)))
    checkTrue(all(genes(gs2) %in% genes(gs12)))
    checkIdentical(genes(GSEABase::setdiff(gs12, gs1)), genes(gs2))
    checkIdentical(genes(GSEABase::setdiff(gs12, gs2)), genes(gs1))
}

test_GCS_LogicalOverlapping <- function() {
    gcss <- .colorBroadSets()
    gs1 <- gcss[[1]]

    idx1 <- sample(seq_along(genes(gs1)), 20)
    idx2 <- sample(seq_along(genes(gcss[[2]])), 20)

    gs2 <-
        GeneColorSet(type=setType(gs1),
                     setName="123", setIdentifier="456",
                     genes=c(genes(gs1)[idx1], genes(gcss[[2]])[idx2]),
                     phenotype=phenotype(gs1),
                     geneColor=factor(c(
                       as.character(geneColor(gs1))[idx1],
                       as.character(geneColor(gcss[[2]]))[idx2])),
                     phenotypeColor=factor(c(
                       as.character(phenotypeColor(gs1))[idx1],
                       as.character(phenotypeColor(gcss[[2]]))[idx2])))
    checkTrue(all(genes(gs1 | gs2) %in%
                  genes(GSEABase::setdiff(gs1, gs2) |
                        (gs1 & gs2) |
                        GSEABase::setdiff(gs2, gs1))))
}
