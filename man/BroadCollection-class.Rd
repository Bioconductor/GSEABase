\name{BroadCollection-class}
\docType{class}
\alias{BroadCollection-class}
\alias{show,BroadCollection-method}

\title{Class "BroadCollection"}
\description{A \code{\link[CollectionType-class]{CollectionType}} used
  to tag gene sets constructed from files containing Broad XML.}
\section{Objects from the Class}{
  Objects can be created with \code{\link{BroadCollection}}; usually this is
  done during a call to \code{\link{GeneSet}},
  \code{\link{GeneColorSet}}, or \code{\link{getBroadSets}}.
}
\section{Slots}{
  \describe{
    \item{\code{category}:}{Object of class \code{"ScalarCharacter"}
      containing terms from the Broad list of categories, or \code{NA}}
    \item{\code{subCategory}:}{Object of class \code{"ScalarCharacter"}
      containing Broad subcategories, or \code{NA}}
    \item{\code{type}:}{Object of class \code{"ScalarCharacter"}
      containing a string identifier of this class.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{CollectionType}}, directly.
}
\section{Methods}{
  No methods defined with class "BroadCollection" in the signature.
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{
  \code{\linkS4class{CollectionType}}

  Broad gene sets are described at \url{http://www.broad.mit.edu/gsea/}.
}
\examples{
example(CollectionType)
example(getBroadSets)
}
\keyword{classes}
