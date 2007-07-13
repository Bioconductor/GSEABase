.CONSTRUCTORS_CollectionType <- 
    names(getSubclasses(getClass("CollectionType")))

.constructors_Simple(.CONSTRUCTORS_CollectionType[.CONSTRUCTORS_CollectionType != "BroadCollection"])

.GETTERS_CollectionType <- c(collectionType="type")

.getters("CollectionType", .GETTERS_CollectionType)

setMethod("show",
          signature=signature(object="CollectionType"),
          function(object) {
              cat("collectionType:", collectionType(object), "\n")
          })

## BroadCollection

BroadCollection <- function(category="c1", subCategory=NA, ...) {
    if (length(category)!=1 ||
        !(category %in% c("c1", "c2", "c3", "c4")))
        stop(sprintf("invalid BroadCollection category: '%s'",
                     paste(category, collapse="', '")))
    new("BroadCollection",
        category=mkScalar(category),
        subCategory=mkScalar(as.character(subCategory)))
}

.GETTERS_BroadCollection <-
    c(bcCategory="category", bcSubCategory="subCategory")

.getters("BroadCollection", .GETTERS_BroadCollection)

setMethod("show",
          signature=signature(object="BroadCollection"),
          function(object) {
              cat("collectionType: ", collectionType(object), "\n",
                  "  bcCategory: ",
                  switch(bcCategory(object),
                         c1="c1 (Positional)",
                         c2="c2 (Curated)",
                         c3="c3 (Motif)",
                         c4="c4 (Computational)"), "\n",
                  "  bcSubCategory:  ", bcSubCategory(object), "\n", sep="")
          })

## GOCollection

GOCollection <- function(goIds=as.character(NA),
                         evidenceCode="ANY") {
    codes <- c("IMP", "IGI", "IPI", "ISS", "IDA", "IEP", "IEA", "TAS",
               "NAS", "ND", "IC", "ANY", NA)
    codeOk <- evidenceCode %in% codes
    if (!all(codeOk))
        stopf("evidenceCode invalid: '%s'",
              paste(evidenceCode[!codeOk], collapse="', '"))
    if ("ANY" %in% evidenceCode)
        evidenceCode <- codes[!codes %in% c("ANY", NA)]
    new("GOCollection",
        goIds=goIds,
        evidenceCode=evidenceCode)
}

.GETTERS_GOCollection <- c("goIds", "evidenceCode")

.getters("GOCollection", .GETTERS_GOCollection)

setMethod("show",
          signature=signature(object="GOCollection"),
          function(object) {
              cat("collectionType:", collectionType(object), "\n",
                  "  goIds:", goIds(object), "\n",
                  "  evidenceCode:", evidenceCode(object), "\n")
          })
