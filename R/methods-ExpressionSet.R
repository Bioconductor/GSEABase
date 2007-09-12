setMethod("[",
          signature(x="ExpressionSet", i="GeneSet"),
          function(x, i, j, ..., drop=TRUE) {
              if (!identical(geneIdType(i),
                             AnnotationIdentifier(annotation(x))))
                i <- mapIdentifiers(i, AnnotationIdentifier(annotation(x)))
              genesIds <- geneIds(i)
              ogenes <- genesIds[genesIds %in% featureNames(x)]
              if (missing(j))
                  x[ogenes,...]
              else
                  x[ogenes, j, ...]
          })
