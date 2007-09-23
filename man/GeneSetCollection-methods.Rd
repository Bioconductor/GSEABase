\name{GeneSetCollection-methods}
\docType{methods}
\alias{GeneSetCollection}
\alias{GeneSetCollection-methods}
\alias{GeneSetCollection,GeneSet,missing,missing-method}
\alias{GeneSetCollection,list,missing,missing-method}
\alias{GeneSetCollection,ExpressionSet,missing,CollectionType-method}
\alias{GeneSetCollection,ExpressionSet,missing,CollectionIdType-method}
\alias{GeneSetCollection,missing,AnnotationIdentifier,CollectionType-method}
\alias{GeneSetCollection,missing,AnnotationIdentifier,CollectionIdType-method}
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
      \code{signature(object = "missing", idType = "AnnotationIdentifier", setType = "CollectionType")},
      \code{signature(object="missing", idType="AnnotationIdentifier", setType="CollectionIdType")}
    }{
      Construct a gene set collection of \code{CollectionType} entities (e.g., pathways for
      \code{KEGGCollection}, protein families for \code{PfamCollection})
      implied by the map found in \code{annotation(idType)}. If \code{setType} is a
      \code{CollectionIdType} and \code{length(ids(setType))>0}, the
      gene set collection is filtered to contain only those sets implied
      by the ids.
    }

    \item{
      \code{signature(object = "ExpressionSet", idType = "missing", setType = "CollectionType")},
      \code{signature(object="ExpressionSet", idType="missing", setType="CollectionIdType")}
    }{
      Construct a gene set collection using the \code{annotation} and
      \code{featureNames} of \code{object} to identify elements for 
      \code{CollectionType} gene sets (e.g., pathways for
      \code{KEGGCollection}, protein families for \code{PfamCollection})
      implied by \code{object}. The gene set
      collection contains only those \code{AnnotationIdentifiers} found
      in \code{featureNames(object)}; if \code{setType} is a
      \code{CollectionIdType} and \code{length(ids(setType))>0}, the
      gene set collection is further filtered to contain only those sets
      implied by the ids.
    }

    \item{
      \code{signature(object="ExpressionSet", idType="missing", setType="GOCollection")}
    }{
      Construct a gene set collection using the \code{annotation} and
      \code{featureNames} of \code{object} to identify 
      \code{GO} pathways implied by \code{object}. The map between
      \code{featureNames} and \code{GO} pathway identifiers is derived
      from the \code{GO2PROBE} table of the annotation
      package of \code{object}. The gene set collection contains only
      those \code{AnnotationIdentifiers} found in
      \code{featureNames(object)}. The \code{evidenceCode} of
      \code{GOCollection} can be used to restrict the pathways seleted to
      those with matching evidence codes.
    }

    \item{
      \code{signature(object="missing", idType="AnnotationIdentifier", setType="GOCollection")}
    }{
      Construct a gene set collection containing all GO pathways
      referenced in the annotation package of \code{idType}. See the
      \code{signature(object="ExpressionSet", idType="missing",
        setType="GOCollection")} method for details.
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
