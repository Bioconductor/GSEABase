\name{GeneSetCollection-class}
\docType{class}
\alias{GeneSetCollection-class}
\alias{[<-,GeneSetCollection,ANY,ANY,ANY-method}
\alias{[<-,GeneSetCollection,ANY,ANY,GeneSet-method}
\alias{[<-,GeneSetCollection,character,ANY,GeneSet-method}
\alias{[,GeneSetCollection,logical-method}
\alias{[,GeneSetCollection,numeric-method}
\alias{[,GeneSetCollection,character-method}
\alias{[[<-,GeneSetCollection,ANY,ANY,ANY-method}
\alias{[[<-,GeneSetCollection,character,ANY,GeneSet-method}
\alias{[[,GeneSetCollection,character-method}
\alias{names,GeneSetCollection-method}

\title{Class "GeneSetCollection"}
\description{
  a \code{GeneSetCollection} is a collection of related
  \code{\linkS4class{GeneSet}}s. The collection can mix and match
  different types of gene sets. Members of the collection are refered to
  by the \code{setName}s of each gene set.
}
  
\section{Objects from the Class}{
  Construct a \code{GeneSetCollection} with a
  \code{\link{GeneSetCollection}} method, e.g., from a list of gene sets
  or with several gene sets provided as argument to the constructor. See
  examples below.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"}, containing the
      gene sets.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\link[Biobase:class.AssayData]{AssayData}"}, by class "list", distance 2.
}
\section{Methods}{
  Gene set collection construction
  \describe{
    \item{GeneSetCollection}{See \code{\link{GeneSetCollection}} methods
      and \code{getBroadSets} for convenient construction methods.}
  }

  Collection access and subsetting (operations on lists, such as
  \code{length}, also work on \code{GeneSetCollection}).

  \describe{
    \item{names}{\code{signature(x = "GeneSetCollection")}: return the
      \code{setName} of each gene set in the colloection.}
    \item{[<-}{
      \code{signature(x = "GeneSetCollection", i = "ANY", j = "ANY",
        value = "ANY")},
      \code{signature(x = "GeneSetCollection", i = "ANY", j = "ANY",
        value = "GeneSet")},
      \code{signature(x = "GeneSetCollection", i = "character", j =
        "ANY", value = "GeneSet")}: assign new sets to existing set
      members. To add entirely new sets, use a
      \code{\link{GeneSetCollection}} constructor.
    }
    \item{[}{
      \code{signature(x = "GeneSetCollection", i = "logical")}, 
      \code{signature(x = "GeneSetCollection", i = "numeric")},
      \code{signature(x = "GeneSetCollection", i = "character")}: create
      a \code{GeneSetCollection} consisting of a subset of the current
      set. All indicies \code{i} must already be present in the set.
    }
    \item{[[<-}{
      \code{signature(x = "GeneSetCollection", i = "ANY", j = "ANY",
        value = "ANY")},
      \code{signature(x = "GeneSetCollection", i = "character", j =
        "ANY", value = "GeneSet")},
      \code{signature(x = "GeneSetCollection", i = "character")}: Select
      a single gene set from the collection.
    }
  }
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{
  \code{\linkS4class{GeneSet}}, \code{\linkS4class{GeneColorSet}}.
}
\examples{
gs1 <- GeneSet(setName="set1", setIdentifier="101")
gs2 <- GeneSet(setName="set2", setIdentifier="102")

## construct from indivdiual elements...
gsc <- GeneSetCollection(gs1, gs2)
## or from a list
gsc <- GeneSetCollection(list(gs1, gs2))

## 'names' are the setNames
names(gsc)

## a collection of a single gene set
gsc["set1"]
## a gene set
gsc[["set1"]]

## set names must be unique
try(GeneSetCollection(gs1, gs1))
try(gsc[c("set1", "set1")])
}
\keyword{classes}
