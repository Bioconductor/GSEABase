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
              cat("collectionType:",
                  collectionType(object),
                  paste("(category: ",
                        switch(bcCategory(object),
                               c1="c1, Positional",
                               c2="c2, Curated",
                               c3="c3, Motif",
                               c4="c4, Computational"),
                        " / ", bcSubCategory(object), ")", sep=""),
                  "\n")
          })
