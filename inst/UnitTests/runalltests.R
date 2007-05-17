library("GSEABase")
library("RUnit")

options(warn=1)

dirs <- '.'

testFilePat <- ".*_test\\.R$"

allSuite <- defineTestSuite(name="allSuite",
                            dirs=dirs,
                            testFileRegexp=testFilePat,
                            rngKind="default",
                            rngNormalKind="default")

testData <- runTestSuite(allSuite)
printTextProtocol(testData, showDetails=FALSE)

