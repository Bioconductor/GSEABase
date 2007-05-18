test_DefaultConstructor <- function() {
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
    checkEquals(new("ScalarCharacter", "unique!"), setIdentifier(gs))
    checkEquals(new("ScalarCharacter", "TestSet"), geneSetName(gs))
}
