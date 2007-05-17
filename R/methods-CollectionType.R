setMethod("collectionType",
          signature=signature(object="CollectionType"),
          function(object) object@type)

setMethod("show",
          signature=signature(object="CollectionType"),
          function(object) {
              cat("collectionType:", collectionType(object), "\n")
          })

setMethod("show",
          signature=signature(object="BroadCollection"),
          function(object) {
              callNextMethod()
              cat("category and sub-category: ",
                  object@category, ", ",
                  object@subCategory, "\n",
                  sep="")
          })
