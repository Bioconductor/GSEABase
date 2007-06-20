\name{GeneColorSet}
\docType{methods}
\alias{GeneColorSet}
\alias{GeneColorSet-methods}
\alias{GeneColorSet,missing-method}
\alias{GeneColorSet,character-method}
\alias{GeneColorSet,GeneIdentifierType-method}
\alias{GeneColorSet,ExpressionSet-method}
\alias{GeneColorSet,GeneSet-method}
\title{Methods to Construct "GeneColorSet" Instances}
\description{
  \code{GeneColorSet} is a generic for constructing gene color sets
  (i.e., gene sets with "coloring" to indicate how features of genes and
  phenotypes are associated).  }
\section{Methods}{

  Available methods are the same as those for
  \code{GeneSet}, but a \code{GeneColorSet} requires an additional
  \code{phenotype} argument to identify the phenotype that is being
  colored. See documentation for \code{\linkS4class{GeneColorSet}} for
  examples.

  An additional method is:
  \describe{
    \item{\code{signature(type = "GeneSet",
        phenotype="character")}}{This method constructs a 'color' gene
      set from an uncolored gene set.}
  }
}
\seealso{
  \code{\link{GeneColorSet-class}}
}
\keyword{methods}
