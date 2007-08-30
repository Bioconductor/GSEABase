\name{CollectionType}
\alias{CollectionType}

\alias{NullCollection}
\alias{ExpressionSetCollection}
\alias{KEGGCollection}
\alias{OMIMCollection}
\alias{PMIDCollection}
\alias{ChrCollection}
\alias{ChrlocCollection}
\alias{MapCollection}
\alias{GOCollection}
\alias{BroadCollection}

\title{Collection Type Class Constructors}
\description{
     These functions construct collection types.
     Collection types can be used in manipulating (e.g., selecting)
     sets, and can contain information specific to particular sets
     (e.g., 'category' and 'subcategory' classifications of
     'BroadCollection'.)
}
\usage{
NullCollection(...)
ExpressionSetCollection(...)
KEGGCollection(ids,...)
OMIMCollection(ids,...)
PMIDCollection(ids,...)
ChrCollection(ids,...)
ChrlocCollection(ids,...)
MapCollection(ids,...)
GOCollection(ids, evidenceCode="ANY", ...)
BroadCollection(category, subCategory=NA, ...)
}

\arguments{
  \item{category}{(Required) Broad category, one of "c1" (postitional), "c2"
    (curated), "c3" (motif), "c4" (computational).}
  \item{subCategory}{(Optional) Sub-category; no controlled vocabulary.}
  \item{ids}{(Optional) Character vector of identifiers (e.g., GO, KEGG,
    or PMID terms).}
  \item{evidenceCode}{(Optional) Character vector of GO evidence codes
    to be included, or "ANY" (any identifier; the default).}
  \item{...}{Additional arguments, usually none but see specific
    \code{linkS4class{CollectionType}} classes  for possibilities.}
}

\value{
  An object of the same class as the function name, initialized as
  appropriate for the collection.
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{
  \code{\linkS4class{CollectionType}}, 
}
\examples{
NullCollection()

## NullCollection when no collection type specified
collectionType(GeneSet())
collectionType(GeneSet(collectionType=GOCollection()))

## fl could be a url
fl <- system.file("extdata", "Broad.xml", package="GSEABase")
gs1 <- getBroadSets(fl)[[1]]
collectionType(gs1) # BroadCollection

## new BroadCollection, with different category
bc <- BroadCollection(category=new("ScalarCharacter","c2"))
## change collectionType of gs2
gs2 <- gs1
collectionType(gs2) <- NullCollection()
}
\keyword{manip}
