test_ES_subset <- function() {
    data(sample.ExpressionSet)
    gss <- getBroadSets(system.file("extdata", "Broad.xml",
                                    package="GSEABase"))
    es <- sample.ExpressionSet[gss[[2]],]
    checkEquals(c(Features=1), nrow(es))
    m <- mapIdentifiers(gss[[2]], AnnotationIdentifier(annotation(es)))
    checkTrue(all(featureNames(es) %in% geneIds(m)))
    checkTrue(all.equal(es, sample.ExpressionSet[m,]))
}
