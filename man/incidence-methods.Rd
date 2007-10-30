\name{incidence-methods}
\docType{methods}
\alias{incidence}
\alias{incidence-methods}
\alias{incidence,GeneSet-method}
\alias{incidence,GeneSetCollection-method}
\title{Methods for Constructing Incidence Matricies Between GeneSets}
\description{
  An incidence matrix summarizes shared membership of gene identifiers
  across (pairs of) gene sets.
}
\section{Methods}{
\describe{
  The return value is a matrix with rows representing gene sets and
  columns genes.

  \item{
    \code{signature(x="GeneSet", ...)}
    \code{signature(x="GeneColorSet", ...)}
  }{
    All additional arguments \code{...} are of the same class as
    \code{x}. The incidence matrix contains elements 0 (genes not
    present) or 1 (genes present).
  }

  \item{
    \code{signature(x = "GeneSetCollection", ...)}
  }{
    Additional arguments \code{...} can be of class
    \code{GeneSetCollection} or \code{GeneSet}. The incidence matrix
    contains elements 0 (genes not present) or 1 (genes present).
  }
}}
\examples{
fl <- system.file("extdata", "Broad.xml", package="GSEABase")
gss <- getBroadSets(fl) # GeneSetCollection of 2 sets
## From one or more GeneSetCollections...
imat <- incidence(gss)
dim(imat)
imat[,c(1:3,ncol(imat)-3+1:3)]

## .. or GeneSets
imat1 <- incidence(gss[[1]], gss[[2]], gss[[1]])
imat1[,1:5]
}
\keyword{methods}
