\name{GeneSetCollection-methods}
\docType{methods}
\alias{GeneSetCollection}
\alias{GeneSetCollection-methods}
\alias{GeneSetCollection,GeneSet-method}
\alias{GeneSetCollection,list-method}
\title{Methods to construct GeneSetCollection instances}
\description{
  Use \code{GeneSetCollection} to construct a collection of gene sets
  from \code{\linkS4class{GeneSet}} arguments, or a list of
  \code{GeneSet}s.
}
\usage{
GeneSetCollection(object, ...)
}
\arguments{
  \item{object}{An argument determining how the gene set collection will
    be created, as described in the methods section.}
  \item{...}{Additional arugments for gene set collection construction,
    as described below.}
}
\section{Methods}{
  \describe{
    \item{object = "GeneSet"}{Construct a gene set collection from one
      or more \code{GeneSet} arugments.}
    \item{object = "list"}{Construct a gene set collection from a list
      of \code{GeneSet}s.}
  }
}
\seealso{
  \code{linkS4class{GeneSetCollection}}
}
\examples{
gs1 <- GeneSet(setName="set1", setIdentifier="101")
gs2 <- GeneSet(setName="set2", setIdentifier="102")

## construct from indivdiual elements...
gsc <- GeneSetCollection(gs1, gs2)
## or from a list
gsc <- GeneSetCollection(list(gs1, gs2))

## set names must be unique
try(GeneSetCollection(gs1, gs1))
}
\keyword{methods}
