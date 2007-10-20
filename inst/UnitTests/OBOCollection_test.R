fl <- system.file("extdata", "goslim_plant.obo", package="GSEABase")
obo <- getOBOCollection(fl)

test_getOBOCollection <- function() {
    obj <- getOBOCollection(fl, evidenceCode="TAS")
    checkEquals(106, length(ids(obj)))
    checkEquals("TAS", evidenceCode(obj))
}

test_obo_subsets <- function() {
    s <- c("Generic GO slim", "GOA and proteome slim",
           "Plant GO slim", "Yeast GO slim", "Prokaryotic GO subset")
    names(s) <- c("goslim_generic", "goslim_goa",
                  "goslim_plant", "goslim_yeast", "gosubset_prok")

    checkIdentical(s, subsets(obo))
    checkIdentical(s, subsets(obo, "named"))
    checkIdentical(names(s), subsets(obo, "key"))
    checkIdentical(as.vector(s), subsets(obo, "value"))
    checkIdentical(paste(names(s), " (", s, ")", sep=""),
                   subsets(obo, "full"))
}

test_obo_subset <- function() {
    checkIdentical(obo, obo[])

    obo0 <- obo[evidenceCode="TAS"]
    checkIdentical("TAS", evidenceCode(obo0))
    checkIdentical(ids(obo), ids(obo0))

    obo1 <- obo["goslim_yeast"]
    checkIdentical("goslim_yeast", subsets(obo1, "key"))
    checkEquals(46, length(ids(obo1)))

    obo2 <- obo[c("goslim_yeast", "gosubset_prok")]
    checkIdentical(c("goslim_yeast", "gosubset_prok"),
                   subsets(obo2, "key"))
    checkEquals(90, length(ids(obo2)))

    obo3 <- obo[c("goslim_yeast", "gosubset_prok"),
                evidenceCode="TAS"]
    checkIdentical(ids(obo2), ids(obo3))
    checkIdentical("TAS", evidenceCode(obo3))

    checkException(obo[ids="GO:0000003"],
                   msg="ids must match those implied by subsets",
                   silent=TRUE)
    checkTrue(validObject(obo[ids=ids(obo)]))
}
