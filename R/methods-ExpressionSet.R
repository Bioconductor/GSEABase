setMethod("[",
          signature(x="ExpressionSet", i="GeneSet"),
          function(x, i, j, ..., drop=TRUE) {
              geneIdType(i) <- AnnotationIdentifier(annotation(x))
              genesIds <- geneIds(i)
              ogenes <- genesIds[genesIds %in% featureNames(x)]
              if (missing(j))
                  x[ogenes,...]
              else
                  x[ogenes, j, ...]
          })
