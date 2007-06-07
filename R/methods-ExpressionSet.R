setMethod("[",
          signature(x="ExpressionSet", i="GeneSet"),
          function(x, i, j, ..., drop=TRUE) {
              setType(i) <- AnnotationIdentifier(annotation(x))
              genes <- genes(i)
              ogenes <- genes[genes %in% featureNames(x)]
              if (missing(j))
                  x[ogenes,...]
              else
                  x[ogenes, j, ...]
          })
