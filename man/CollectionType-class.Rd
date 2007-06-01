\name{CollectionType-class}
\docType{class}
% CollectionType-class
\alias{CollectionType-class}
\alias{collectionType,CollectionType-method}
\alias{show,CollectionType-method}
% AdHocCollection-class
\alias{AdHocCollection-class}
% ExpressionSetCollection-class
\alias{ExpressionSetCollection-class}
% BroadCollection-class
\alias{BroadCollection-class}
\alias{bcCategory,BroadCollection-method}
\alias{bcSubCategory,BroadCollection-method}
\alias{show,BroadCollection-method}
% GOCollection-class
\alias{GOCollection-class}
% KEGGCollection-class
\alias{KEGGCollection-class}

\title{Class "CollectionType"}

\description{
  These classes provides a way to tag the origin of a
  \code{GeneSet}. Collection types can be used in manipulating
  (e.g., selecting) sets, and can contain information specific to
  particular sets (e.g., \code{category} and \code{subcategory}
  classifications of \code{BroadCollection}.)
}

\section{Objects from the Class}{
  The following classes can tag gene sets:
  \describe{
    \item{AdHocCollection}{No formal collection information available.}
    \item{BroadCollection}{Derived from, or destined to, Broad
      XML. Usually created and written \code{\link{getBroadSets}},
      \code{\link{toBroadXML}}.}
    \item{ExpressionSetCollection}{Derived from
      \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. Usually
      created during a call to \code{\link{GeneSet}} or
      \code{\link{GeneColorSet}}.}
    \item{GOCollection}{Collection derived using Gene Ontology (GO)
      terms.}
    \item{KEGGCollection}{Collection derived using KEGG terms.}
  }

  Objects are instantiated with calls to \code{\link{CollectionType}}
  constructors.
}

\section{Slots}{
  \code{CollectionType} classes have the slot:
  \describe{
    \item{\code{type}:}{Object of class \code{"ScalarCharacter"}
      containing the character string representation of this
      \code{CollectionType}.}
  }

  \code{BroadCollection} has slots:
  \describe{
     \item{\code{category}:}{Object of class \code{"ScalarCharacter"}
      containing terms from the Broad list of categories, or \code{NA}}
    \item{\code{subCategory}:}{Object of class \code{"ScalarCharacter"}
      containing Broad sub-categories, or \code{NA}}
  }
}

\section{Methods}{
  \code{CollectionType} classes have methods:
  \describe{
    \item{collectionType<-}{\code{signature(object = "GeneSet", value =
        "CollectionType")}: Replace the \code{CollectionType} }
    \item{collectionType}{\code{signature(object = "CollectionType")}:
      Retrieve the collection type.}
    \item{show}{\code{signature(object = "CollectionType")}: display the
      collection type.}
  }

  \code{BroadCollection} has methods:
  \describe{
    \item{bcCategory}{Retrieve the \code{category} of the Broad
      collection.}
    \item{bcSubCategory}{Retrieve the sub-category of the Broad
      collection.}
  }
}

\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{
  \code{\link{CollectionType}} consturctors.
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
