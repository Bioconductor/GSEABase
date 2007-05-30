\name{CollectionType-class}
\docType{class}

\alias{CollectionType-class}

\alias{collectionType,CollectionType-method}
\alias{show,CollectionType-method}

\title{Class "CollectionType"}
\description{
  This class provides a way to tag the origin of a \code{GeneSet}. For
  instance, a \code{GeneSet} derived from an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}
  might have \code{CollectionType} of \code{ExpressionSetCollection}, whereas
  a gene set derived from Broad XML has
  \code{BroadCollection}. Collection types can be used in mannipulating
  (e.g., selecting) sets, and can contain information specific to
  particular sets (e.g., \code{category} and \code{subcategory}
  classifications of \code{BroadCollection}.)} 
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
	 \describe{
    \item{\code{type}:}{Object of class \code{"ScalarCharacter"}
      containing the character string representation of this
      \code{CollectionType}.}
  }
}
\section{Methods}{
  \describe{
    \item{collectionType<-}{\code{signature(object = "GeneSet", value =
        "CollectionType")}: Replace the \code{CollectionType} }
    \item{collectionType}{\code{signature(object = "CollectionType")}:
      Retrieve the collection type.}
    \item{show}{\code{signature(object = "CollectionType")}: display the
      collection type.}
  }
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{
  The example below lists \code{CollectionType} classes defined in this
  package. See the help pages of these classes for specific information.

}
\examples{
getSubclasses(getClass("CollectionType"))

## Create a CollectionType and ask for its type
collectionType(ExpressionSetCollection())

\dontrun{

## Read two GeneSets from a Broad XML file into a list, verify that
## they are both BroadCollection's. Category / subcategory information
## is unique to Broad collections.
fl <- system.file("extdata", "Broad.xml", package="GSEABase")
sets <- getBroadSets(fl)
sapply(sets, collectionType)
}

## ExpressionSets are tagged with ExpressionSetCollection; there is no
## 'category' information.
data(sample.ExpressionSet)
gs <- GeneSet(sample.ExpressionSet[100:109],
              setName="sample.GeneSet", setIdentifier="123")
collectionType(gs)
}
\keyword{classes}
