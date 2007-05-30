\name{Identifiers}
\alias{<name>Identifier}
\alias{NullIdentifier}
\alias{AnnotationIdentifier}
\alias{EntrezIdentifier}
\alias{SymbolIdentifier}
\alias{PfamIdentifier}
\title{Gene Identifier Class Constructors}
\description{
  These functions construct gene identifiers. A gene identifier is used
  by \code{\linkS4class{GeneSet}} classes to indicate what the list of
  genes represents (e.g., Entrez gene identifiers, Bioconductor
  annotations).
}
\usage{
NullIdentifier(...)
AnnotationIdentifier(...)
EntrezIdentifier(...)
SymbolIdentifier(...)
PfamIdentifier(...)
}
\arguments{
  \item{\dots}{Additional arguments, usually none but see specific
    identifier pages for possibilities.}
}
\value{
  An object of the same class as the function name, initialized as
  appropriate for the identifier.
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{
  \code{\linkS4class{NullIdentifier}}, 
  \code{\linkS4class{AnnotationIdentifier}}, 
  \code{\linkS4class{EntrezIdentifier}}, 
  \code{\linkS4class{SymbolIdentifier}}, 
  \code{\linkS4class{PfamIdentifier}}
}
\examples{
NullIdentifier()

data(sample.ExpressionSet)
gs1 <- GeneSet(sample.ExpressionSet[100:109],
               setName="sample1", setIdentifier="100")
setType(gs1) # AnnotationIdentifier

genes <- featureNames(sample.ExpressionSet)[100:109]
gs2 <- GeneSet(genes=genes,
               setName="sample1", setIdentifier="101")
setType(gs2) # NullIdentifier, since no info about genes provided

## Convert...
setType(gs2) <- AnnotationIdentifier()
setType(gs2)
## ...or provide more explicit construction
gs3 <- GeneSet(genes=genes, type=AnnotationIdentifier(),
               setName="sample1", setIdentifier="102")
}
\keyword{manip}
