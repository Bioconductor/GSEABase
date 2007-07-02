\name{GeneSet}
\docType{methods}
\alias{GeneSet}
\alias{GeneSet-methods}
\alias{GeneSet,missing-method}
\alias{GeneSet,character-method}
\alias{GeneSet,GeneIdentifierType-method}
\alias{GeneSet,ExpressionSet-method}
\title{Methods to construct GeneSet instances}
\description{
  Use \code{GeneSet} to construct gene sets from \code{ExpressionSet},
  character vector, or other objects.
}
\usage{
GeneSet(type, ...)
}
\arguments{
  \item{type}{An argument determining how the gene set will be created,
    as described in the Methods section.}
  \item{...}{Additional arguments for gene set construction. 
    Methods have required arguments, as outlined below; additional
    arguments correspond to slot names \code{\linkS4class{GeneSet}}.}
}
\section{Methods}{
  \describe{
    \item{
      \code{signature(type = "missing", ...)}}{
      Construct an empty gene set.}

    \item{
      \code{signature(type = "character", ...)}}{
      Construct a gene set using identifiers \code{type}.}

    \item{
      \code{signature=(type = "GeneIdentifierType", ...)}}{
      Construct an empty gene set. The gene set has \code{geneIdType}
      created from the \code{GeneIdentifierType} of \code{type}.}

    \item{
      \code{signature=(type = "ExpressionSet", ...)}}{
      Construct a gene set from an
      \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. \code{geneIdType}
      is set to \code{AnnotationIdentifier}; the annotation field and
      annotation package of the \code{ExpressionSet} are consulted to
      determine \code{organism}, if possible. Short and long
      descriptions from the \code{ExpressionSet} \code{experimentData}
      title and abstract; pub med ids, urls, and contributor are also
      derived from \code{experimentData}.}
  }
}
\seealso{
  \code{\link{GeneSet-class}}
  \code{\link{GeneColorSet-class}}
}
\examples{
## Empty gene set
GeneSet()

## Gene set from ExpressionSet
data(sample.ExpressionSet)
gs1 <- GeneSet(sample.ExpressionSet[100:109])

## GeneSet from Broad XML; 'fl' could be a url
fl <- system.file("extdata", "Broad.xml", package="GSEABase")
gs2 <- getBroadSets(fl)[[1]] # actually, a list of two gene sets

## GeneSet from list of gene identifiers
geneIds <- geneIds(gs2) # any character vector would do
gs3 <- GeneSet(geneIds)
## unspecified set type, so...
is(geneIdType(gs3), "NullIdentifier") == TRUE
## update set type to match encoding of identifiers
geneIdType(gs2)
geneIdType(gs3) <- EntrezIdentifier()
## other ways of accomplishing the same
gs4 <- GeneSet(geneIds, geneIdType=EntrezIdentifier())
gs5 <- GeneSet(EntrezIdentifier(), geneIds=geneIds)
}
\keyword{methods}
