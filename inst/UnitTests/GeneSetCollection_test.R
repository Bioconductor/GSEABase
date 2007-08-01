.broadSets <- function()
    getBroadSets(system.file("extdata", "Broad.xml", package="GSEABase"))

.gsc <- function() {
    gs1 <- GeneSet(setName="set1", setIdentifier="id1",
                   geneIds=LETTERS[1:5])
    gs2 <- GeneSet(setName="set2", setIdentifier="id2",
                   geneIds=letters[1:5])
    GeneSetCollection(list(gs1, gs2))
}

test_GSC_list_constructor <- function() {
    gs1 <- GeneSet(setName="123", setIdentifier="456",
                   geneIds=LETTERS[1:5])
    gs2 <- GeneSet(setName="234", setIdentifier="567",
                   geneIds=letters[1:5])
    gsc <- GeneSetCollection(list(gs1, gs2))
    checkTrue(validObject(gsc))
    checkEquals(2, length(gsc))
    checkIdentical(gs1, gsc[[1]])
    checkIdentical(gs2, gsc[[2]])

    ## duplicate entries
    checkException(GeneSetCollection(gs1, gs1), silent=TRUE)
}

test_GSC_docs_constructor <- function() {
    gs1 <- GeneSet(setName="123", setIdentifier="456",
                   geneIds=LETTERS[1:5])
    gs2 <- GeneSet(setName="234", setIdentifier="567",
                   geneIds=letters[1:5])
    gsc <- GeneSetCollection(gs1, gs2)
    checkTrue(validObject(gsc))
    checkEquals(2, length(gsc))
    checkIdentical(gs1, gsc[[1]])
    checkIdentical(gs2, gsc[[2]])
}

test_GSC_validity <- function() {
    gsc <- .gsc()
    gsc@.Data <- append(gsc@.Data, 1)
    checkException(validObject(gsc), silent=TRUE)
}

test_GSC_length <- function() {
    checkTrue(length(.gsc())==2)
}

test_GSC_names <- function() {
    checkTrue(all(c("set1", "set2")==names(.gsc())))
}

test_GSC_subset_by_name<- function() {
    gsc <- .gsc()

    gsc1 <- gsc["set1"]
    checkTrue(is(gsc1, "GeneSetCollection"))
    checkTrue(validObject(gsc1))
    checkEquals(1, length(gsc1))
    checkEquals("set1", names(gsc1))

    gsc21 <- gsc[c("set2", "set1")]
    checkTrue(validObject(gsc21))
    checkEquals(2, length(gsc21))
    checkTrue(all(c("set2", "set1")==names(gsc21)))

    checkException(gsc["set3"], silent=TRUE) # no element
    checkException(gsc[c("set1", "set1")], silent=TRUE) # duplicate entries
}

test_GSC_subset_by_index <- function() {
    gsc <- .gsc()

    gsc1 <- gsc[1]
    checkTrue(is(gsc1, "GeneSetCollection"))
    checkTrue(validObject(gsc1))
    checkEquals(1, length(gsc1))
    checkEquals("set1", names(gsc1))

    gsc21 <- gsc[2:1]
    checkTrue(validObject(gsc21))
    checkEquals(2, length(gsc21))
    checkTrue(all(c("set2", "set1")==names(gsc21)))

    checkException(gsc[3], silent=TRUE) # no element
    checkException(gsc[c(1,1)], silent=TRUE) # duplicate entries
}

test_GSC_subset_by_logical <- function() {
    gsc <- .gsc()
    checkException(gsc[rep(TRUE, 3)], silent=TRUE) # out-of-bounds
}

test_GSC_subset2 <- function() {
    gsc <- .gsc()
    gsc2 <- gsc[[2]]
    checkTrue(is(gsc2, "GeneSet"))
    checkTrue(validObject(gsc2))
    checkTrue("set2"==setName(gsc2))

    gsc2 <- gsc[["set2"]]
    checkTrue(is(gsc2, "GeneSet"))
    checkTrue(validObject(gsc2))
    checkTrue("set2"==setName(gsc2))

    ## subscript out of bounds
    checkException(gsc[[c(1,2)]], silent=TRUE)
    checkException(gsc[[c("set1", "set2")]], silent=TRUE)
    checkException(gsc[[3]], silent=TRUE)
    checkException(gsc[["set3"]], silent=TRUE)
}

## test_GSC_subset_assign <- function() {
##     checkTrue(FALSE)
## }

## test_GSC_subset2_assign <- function() {
##     checkTrue(FALSE)
## }

test_GSC_incidence <- function() {
    gss <- .broadSets()
    res <- incidence(gss)
    checkTrue(all(dim(res)==c(2, 215)))
    checkTrue(sum(res)== 215)
    res <- incidence(gss, gss)
    checkTrue(all(dim(res)==c(4, 215)))
    checkTrue(sum(res)== 430)
}
