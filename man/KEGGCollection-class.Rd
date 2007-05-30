\name{KEGGCollection-class}
\docType{class}
\alias{KEGGCollection-class}

\title{Class "KEGGCollection"}
\description{A \code{\link[CollectionType-class]{CollectionType}} used
  to tag gene sets constructed from KEGG terms.}
\section{Objects from the Class}{
  Objects can be created with \code{\link{KEGGCollection}}; usually this is
  done during a call to \code{\link{GeneSet}},
  \code{\link{GeneColorSet}}.
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
