\name{GeneSetCollection-class}
\docType{class}
\alias{GeneSetCollection-class}

\alias{updateObject,GeneSetCollection-method}
\alias{names,GeneSetCollection-method}
\alias{geneIds,GeneSetCollection-method}
\alias{geneIds<-,GeneSetCollection,list-method}

\alias{Logic,character,GeneSetCollection-method}
\alias{Logic,GeneSet,GeneSetCollection-method}
\alias{|,GeneSetCollection,GeneSetCollection-method}
\alias{|,GeneSetCollection,ANY-method}
\alias{|,ANY,GeneSetCollection-method}
\alias{&,GeneSetCollection,GeneSetCollection-method}
\alias{&,GeneSetCollection,ANY-method}
\alias{&,ANY,GeneSetCollection-method}
\alias{intersect,ANY,GeneSetCollection-method}
\alias{intersect,GeneSetCollection,ANY-method}
\alias{setdiff,ANY,GeneSetCollection-method}
\alias{setdiff,GeneSetCollection,ANY-method}
\alias{toGmt,GeneSetCollection-method}
\alias{union,ANY,GeneSetCollection-method}
\alias{union,GeneSetCollection,ANY-method}

\alias{[<-,GeneSetCollection,ANY,ANY,ANY-method}
\alias{[<-,GeneSetCollection,ANY,ANY,GeneSet-method}
\alias{[<-,GeneSetCollection,character,ANY,GeneSet-method}
\alias{[,GeneSetCollection,logical,ANY-method}
\alias{[,GeneSetCollection,numeric,ANY-method}
\alias{[,GeneSetCollection,character,ANY-method}
\alias{[[,GeneSetCollection,character-method}
\alias{[[<-,GeneSetCollection,ANY,ANY,ANY-method}
\alias{[[<-,GeneSetCollection,numeric,ANY,GeneSet-method}
\alias{[[<-,GeneSetCollection,character,ANY,GeneSet-method}

\alias{show,GeneSetCollection-method}

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

  Collection access (operations on lists, such as \code{length}, ,
  \code{lapply} also work on \code{GeneSetCollection}).
  \describe{
    \item{geneIds}{\code{signature(object = "GeneSetCollection")}:
      return a list, with each member a character vector of gene
      identifiers from the gene set collection.}
    \item{geneIds<-}{\code{signature(object="GeneSetCollection",
        value="list")}: assign character vectors in \code{value} to
      corresponding \code{geneIds} of \code{object}.}
    \item{names}{\code{signature(x = "GeneSetCollection")}: return the
      \code{setName} of each gene set in the colloection.}
  }

  Logical and subsetting operations
  \describe{
    \item{union}{
      \code{signature(x = "GeneSetCollection", y = "ANY")},
      \code{signature(x = "ANY", y = "GeneSetCollection")}: ...}
    \item{|}{
      \code{signature(e1 = "GeneSetCollection", e2 = "ANY")},
      \code{signautre(e1 = "GeneSet", e2 = "GeneSetCollection")},
      \code{signautre(e1 = "character", e2 = "GeneSetCollection")},
      \code{signature(e1 = "ANY", e2 = "GeneSetCollection")}:
      calculate the logical `or` (union) of all gene identifiers
      in an object over all members of the gene set collection.} 
    \item{intersect}{
      \code{signature(x = "GeneSetCollection", y = "ANY")},
      \code{signature(x = "ANY", y = "GeneSetCollection")}: ...}
    \item{&}{
      \code{signature(e1 = "GeneSetCollection", e2 = "ANY")},
      \code{signautre(e1 = "character", e2 = "GeneSetCollection")},
      \code{signautre(e1 = "GeneSet", e2 = "GeneSetCollection")},
      \code{signature(e1 = "ANY", e2 = "GeneSetCollection")}:
      calculate the logical `and' (intersection) of all gene identifiers
      in a gene set or character vector, over all members of the gene
      set collection.} 
    \item{setdiff}{
      \code{signature(x = "GeneSetCollection", y = "ANY")}:
      calculate the logical set difference betwen all gene sets in a
      collection and the gene identifiers of a gene set or character
      vector. A \code{setdiff} method must be available for
      \code{x="GeneSet"} and the type of \code{y}. }

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
    \item{[[}{
      \code{signature(x = "GeneSetCollection", i = "character")}: 
      Select a single gene set from the collection. Methods for
      \code{i="numeric"} are inherited from \code{list}.}

    \item{[[<-}{
      \code{signature(x = "GeneSetCollection", i = "ANY", j = "ANY", value = "ANY")},
      \code{signature(x = "GeneSetCollection", i = "numeric", j = "ANY", value = "GeneSet")},
      \code{signature(x = "GeneSetCollection", i = "character", j = "ANY", value = "GeneSet")}:
      Replace a gene set in the collecton with another.
      \code{value = "ANY"} serves to stop invalid assignments.} 
  }

  Additional useful methods.

  \describe{
    \item{updateObject}{
      Objects created in previous versions of GSEABase may be
      incompatible with current object definitions. Usually this is
      singalled by an error suggesting that a slot is missing, and a
      recommnedation to use \code{updateObject}. Use \code{updateObject}
      to update a \code{GeneSetCollection} and all contained
      \code{GeneSets} to their current defintion.}

    \item{mapIdentifiers}{
      Convert genes from one \code{GeneIdentifierType} to another. See
      \code{\link{mapIdentifiers}} and specific methods in
      \code{\linkS4class{GeneIdentifierType}} for additional detail.}

    \item{incidence}{Summarize shared membership in genes across gene
      sets. See \code{\link{incidence-methods}}.}

    \item{toGmt}{Export to 'GMT' format file. See \code{\link{toGmt}}.}

    \item{show}{\code{signature(object="GeneSetCollection")}: provide a
      compact representation of \code{object}.}
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
