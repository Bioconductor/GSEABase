\name{details-methods}
\docType{methods}
\alias{details}
\alias{details-methods}
\alias{details,GeneSet-method}
\alias{details,GeneColorSet-method}
\title{Methods for Displaying Detailed GeneSet Information}
\description{
  This generic and methods supplement \code{show}, providing more detail
  on object contents.
}
\section{Methods}{
  \describe{
    Defined methods include:
    \item{
      \code{signature(object = "GeneSet")},
      \code{signature(object = "GeneColorSet")}
    }{
      These methods display information about \code{setIdentifier},
      \code{description}, \code{organism}, \code{pubMedIds},
      \code{urls}, \code{contributor}, \code{setVersion}, and
      \code{creationDate}.
    }
  }
}
\keyword{methods}
