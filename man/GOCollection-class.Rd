\name{GOCollection-class}
\docType{class}
\alias{GOCollection-class}

\title{Class "GOCollection"}
\description{A \code{\linkS4class{CollectionType}} used
  to tag gene sets constructed from GO terms.}
\section{Objects from the Class}{
  Objects can be created with \code{\link{GOCollection}}; usually this is
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
  No methods defined with class "KEGGCollection" in the signature.
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{
  \code{\linkS4class{CollectionType}}
}
\examples{
example(CollectionType)
}
\keyword{classes}
