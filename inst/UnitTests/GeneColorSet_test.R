test_ConstructorNoColorSetArgs <- function() {
    checkException(GeneColorSet(setIdentifier="123",
                                setName="Set name"),
                   silent=TRUE)
}

test_ConstructorDefaultArgs <- function() {
    gs <- GeneColorSet(setIdentifier="123",
                       setName="Set name",
                       phenotype="A phenotype")
    checkTrue(validObject(gs, complete=TRUE))
    checkTrue(length(genes(gs))==0)
    checkTrue(length(geneColor(gs))==0)
    checkTrue(length(phenotypeColor(gs))==0)
}

test_ConstructorAllColorArgs <- function() {
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
    gfactor <-     factor(rep(c("high", "low"), 12))
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
