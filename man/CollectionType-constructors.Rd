\name{Collections}
\alias{<name>Collection}
\alias{AdHocCollection}
\alias{ExpressionSetCollection}
\alias{KEGGCollection}
\alias{GOCollection}
\alias{BroadCollection}
\title{Collection Type Class Constructors}
\description{
  These functions construct gene identifiers. A gene identifier is used
  by \code{\linkS4class{GeneSet}} classes to indicate what the list of
  genes represents (e.g., Entrez gene identifiers, Bioconductor
  annotations).
}
\usage{
AdHocCollection(...)
ExpressionSetCollection(...)
KEGGCollection(...)
GOCollection(...)
BroadCollection(...)
}
\arguments{
  \item{\dots}{Additional arguments, usually none but see specific
    collection pages for possibilities.}
}
\value{
  An object of the same class as the function name, initialized as
  appropriate for the collection.
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{
  \code{\linkS4class{AdHocCollection}}, 
  \code{\linkS4class{ExpressionSetCollection}}, 
  \code{\linkS4class{KEGGCollection}}, 
  \code{\linkS4class{GOCollection}}, 
  \code{\linkS4class{BroadCollection}}
}
\examples{
AdHocCollection()

## Ad hoc when no collection type specified
collectionType(GeneSet(setName="empty", setIdentifier="0"))
collectionType(GeneSet(collectionType=GOCollection(),
                       setName="empty",
                       setIdentifier="0"))

## fl could be a url
fl <- system.file("extdata", "Broad.xml", package="GSEABase")
gs1 <- getBroadSets(fl)[[1]]
collectionType(gs1) # BroadCollection

## new BroadCollection, with different category
bc <- BroadCollection(category=new("ScalarCharacter","c2"))
## change collectionType of gs2
gs2 <- gs1
collectionType(gs2) <- AdHocCollection()
}
\keyword{manip}
