\name{GeneSetCollection-methods}
\docType{methods}
\alias{GeneSetCollection}
\alias{GeneSetCollection-methods}
\alias{GeneSetCollection,GeneSet,missing,missing-method}
\alias{GeneSetCollection,list,missing,missing-method}
\alias{GeneSetCollection,ExpressionSet,missing,KEGGCollection-method}
\alias{GeneSetCollection,ExpressionSet,missing,GOCollection-method}
\alias{GeneSetCollection,missing,AnnotationIdentifier,KEGGCollection-method}
\alias{GeneSetCollection,missing,AnnotationIdentifier,GOCollection-method}
\title{Methods to construct GeneSetCollection instances}
\description{
  Use \code{GeneSetCollection} to construct a collection of gene sets
  from \code{\linkS4class{GeneSet}} arguments, or a list of
  \code{GeneSet}s.
}
\usage{
GeneSetCollection(object, ..., idType, setType)
}
\arguments{
  \item{object}{An argument determining how the gene set collection will
    be created, as described in the methods section.}
  \item{...}{Additional arugments for gene set collection construction,
    as described below.}
  \item{idType}{An argument of class
    \code{\linkS4class{GeneIdentifierType}}, used to
    indicate how the \code{geneIds} will be represented.}
  \item{setType}{An argument of class
    \code{\linkS4class{CollectionType}}, used to indicate how the
    collection is created.}
}
\section{Methods}{
  \describe{
    \item{
      \code{signature(object = "GeneSet",idType="missing", setType="missing")}
    }{
      Construct a gene set collection from one or more \code{GeneSet}
      arugments.
    }
    
    \item{
      \code{signature(object = "list", idType="missing", setType="missing")}
    }{
      Construct a gene set collection from a list of \code{GeneSet}s.
    }

    \item{
      \code{signature(object = "ExpressionSet", idType = "missing", setType = "KEGGCollection")},
      \code{signature(object="ExpressionSet", idType="missing", setType="GOCollection")}
    }{
      Construct a gene set collection using the \code{annotation} and
      \code{featureNames} of \code{object} to identify \code{KEGG} or
      \code{GO} pathways implied by \code{object}. The gene set
      collection contains only those \code{AnnotationIdentifiers} found
      in \code{featurNames(object)}.
    }

    \item{
      \code{signature(object="missing", idType="AnnotationIdentifier", setType="KEGGCollection")},
      \code{signature(object="missing", idType="AnnotationIdentifier", setType="GoCollection")}
    }{
      Construct a gene set collection containing all GO or KEGG pathways
      referenced in the annotation of \code{idType}. The
      \code{evidenceCode} of \code{GOCollection} can be used to restrict
      the pathways seleted to those with matching evidence codes.
    }
  }
}
\seealso{
  \code{\linkS4class{GeneSetCollection}}
}
\examples{
gs1 <- GeneSet(setName="set1", setIdentifier="101")
gs2 <- GeneSet(setName="set2", setIdentifier="102")

## construct from indivdiual elements...
gsc <- GeneSetCollection(gs1, gs2)
## or from a list
gsc <- GeneSetCollection(list(gs1, gs2))

## set names must be unique
try(GeneSetCollection(gs1, gs1))

data(sample.ExpressionSet)
gsc <- GeneSetCollection(sample.ExpressionSet[200:250],
                         setType = GOCollection())
}
\keyword{methods}
