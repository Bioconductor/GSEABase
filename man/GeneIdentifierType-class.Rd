\name{GeneIdentifierType-class}
\docType{class}
\alias{GeneIdentifierType-class}
\alias{setType,GeneIdentifierType-method}
\alias{show,GeneIdentifierType-method}
\alias{validIdentifiers,GeneIdentifierType-method}

\title{Class "GeneIdentifierType"}

\description{This class provides a way to tag the meaning of gene
  symbols in a \code{GeneSet}. For instance, a \code{GeneSet} with gene
  names derived from a Bioconductor \code{annotation} package (e.g., via
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}) initially have a
  \code{GeneIdentifierType} of \code{AnnotationType}.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
	 \describe{
    \item{\code{type}:}{Object of class \code{"ScalarCharacter"}
      containing the character string representation of this
      \code{GeneIdentifierType}.}
  }
}
\section{Methods}{
  \describe{
    \item{GeneSet}{\code{signature(type = "GeneIdentifierType")}:
      Create a new \code{\link{GeneSet}} using identifiers of
      \code{GeneIdentifierType}.}

    \item{GeneColorSet}{\code{signature(type = "GeneIdentifierType")}:
      Create a new \code{\link{GeneColorSet}} using identifiers of
      \code{GeneIdentifierType}.}

    \item{setType<-}{\code{signature(object = "GeneSet", value =
        "GeneIdentifierType")}: Change the \code{GeneIdentifierType} of
      \code{object} to \code{value}, attempting to convert symbols in
      the process.}

    \item{setType}{\code{signature(object = "GeneIdentifierType")}:
      return the type of this \code{object}.}

    \item{show}{\code{signature(object = "GeneIdentifierType")}: display
      this object.}

    \item{validIdentifiers}{\code{signature(identifier =
        "GeneIdentifierType")}: validate that the identiifier are
      legitimate; subclasses are responsible for implementing this
      method.}
  }
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{
   The example below lists \code{GeneIdentifierType} classes defined in
   this package; See the help pages of these classes for specific information.
}
\examples{
getSubclasses(getClass("GeneIdentifierType"))

# create an AnnotationIdentifier, and ask it's type
setType(AnnotationIdentifier())

# Construct a GeneSet from an ExpressionSet, using the 'annotation'
# field of ExpressionSet to recognize the genes as AnnotationType
data(sample.ExpressionSet)
gs <- GeneSet(sample.ExpressionSet[100:109],
              setName="sample.GeneSet", setIdentifier="123")
setType(gs)

## Read a Broad set from the sytem (or a url), and discover their
## GeneIdentifierType
fl <- system.file("extdata", "Broad.xml", package="GSEABase")
bsets <- getBroadSets(fl)
sapply(bsets, setType)

## try to combine gene sets with different set types
try(gs & sets[[1]])

\dontrun{
## Use the annotation package associated with the original
## ExpressionSet to map to EntrezIdentifier() ...
setType(gs) <- EntrezIdentifier()
...## and try again
gs & bsets[[1]]
}

## Create a new identifier
setClass("FooIdentifier",
         contains="GeneIdentifierType",
         prototype=prototype(
           type=new("ScalarCharacter", "Foo")))
## Create a constructor (optional)
FooIdentifier <- function() new("FooIdentifier")
setType(FooIdentifier())

## tidy up
removeClass("FooIdentifier")
}
\keyword{classes}
