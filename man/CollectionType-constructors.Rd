\name{CollectionType}
\alias{CollectionType}

\alias{NullCollection}
\alias{ExpressionSetCollection}
\alias{ComputedCollection}
\alias{KEGGCollection}
\alias{OMIMCollection}
\alias{PfamCollection}
\alias{PrositeCollection}
\alias{PMIDCollection}
\alias{ChrCollection}
\alias{ChrlocCollection}
\alias{MapCollection}
\alias{GOCollection}
\alias{OBOCollection}
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
ComputedCollection(...)
ExpressionSetCollection(...)
ChrCollection(ids,...)
ChrlocCollection(ids,...)
KEGGCollection(ids,...)
MapCollection(ids,...)
OMIMCollection(ids,...)
PMIDCollection(ids,...)
PfamCollection(ids, ...)
PrositeCollection(ids, ...)
GOCollection(ids=character(0), evidenceCode="ANY", ontology="ANY", ..., err=FALSE)
OBOCollection(ids, evidenceCode="ANY", ontology="ANY", ...)
BroadCollection(category, subCategory=NA, ...)
}

\arguments{
  \item{category}{(Required) Broad category, one of "c1" (postitional),
    "c2" (curated), "c3" (motif), "c4" (computational), "c5" (GO), "c6"
    (Oncogenic Pathway Activation Modules) "c7" (Immunologic
    Signatures), "h" (Hallmark).}
  \item{subCategory}{(Optional) Sub-category; no controlled vocabulary.}
  \item{ids}{(Optional) Character vector of identifiers (e.g., GO, KEGG,
    or PMID terms).}
  \item{evidenceCode}{(Optional) Character vector of GO evidence codes
    to be included, or "ANY" (any identifier; the default). 
    Evidence is a property of particular genes, rather than of the
    ontology, so evidenceCode is a convenient way of specifying how
    users of a GOCollection might restrict derived objects (as in
    done during create of a gene set from an expression set).}
  \item{ontology}{(Optional) Character vector of GO ontology terms to be
    included, or "ANY" (any identifier; the default). Unlike evidence
    code, ontology membership is enforced when GOCollection gene sets
    are constructed.}
  \item{err}{(Optional) logical scalar indicating whether non-existent
    GO terms signal an error (\code{TRUE}), or are silently ignored (\code{FALSE}).}
  \item{...}{Additional arguments, usually none but see specific
    \code{\linkS4class{CollectionType}} classes  for possibilities.}
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
bc <- BroadCollection(category="c2")
## change collectionType of gs2
gs2 <- gs1
collectionType(gs2) <- NullCollection()

## OBOCollection
fl <- system.file("extdata", "goslim_plant.obo", package="GSEABase")
getOBOCollection(fl, evidenceCode="TAS") # returns OBOCollection
OBOCollection(c("GO:0008967", "GO:0015119", "GO:0030372", "GO:0002732",
                "GO:0048090"))
}
\keyword{manip}
