.CONSTRUCTORS_CollectionType <- 
    names(getSubclasses(getClass("CollectionType")))

.constructors_Simple(.CONSTRUCTORS_CollectionType)

setMethod("collectionType",
          signature=signature(object="CollectionType"),
          function(object) object@type)

setMethod("show",
          signature=signature(object="CollectionType"),
          function(object) {
              cat("collectionType:", collectionType(object), "\n")
          })

## BroadCollection

setMethod("show",
          signature=signature(object="BroadCollection"),
          function(object) {
              cat("collectionType:",
                  collectionType(object),
                  paste("(category: ", object@category, " / ",
                        object@subCategory, ")", sep=""),
                  "\n")
          })
