\name{ExpressionSetCollection-class}
\docType{class}
\alias{ExpressionSetCollection-class}

\title{Class "ExpressionSetCollection"}
\description{A \code{\link[CollectionType-class]{CollectionType}} used
  to tag gene sets constructed from an ExpressionSet.}
\section{Objects from the Class}{
  Objects can be created with \code{\link{ExpressionSetCollection}}; usually this is
  done during a call to \code{\link{GeneSet}} or \code{\link{GeneColorSet}}.
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"ScalarCharacter"}
      containing a string identifier of this class.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{CollectionType}}, directly.
}
\section{Methods}{
  No methods defined with class "ExpressionSetCollection" in the signature.
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{
  \code{\linkS4class{CollectionType}}
}
\examples{
example(CollectionType)
}
\keyword{classes}
