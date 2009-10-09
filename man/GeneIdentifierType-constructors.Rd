\name{GeneIdentifierType}
\alias{GeneIdentifierType}

\alias{NullIdentifier}
\alias{EnzymeIdentifier}
\alias{ENSEMBLIdentifier}
\alias{GenenameIdentifier}
\alias{RefseqIdentifier}
\alias{SymbolIdentifier}
\alias{UnigeneIdentifier}
\alias{EntrezIdentifier}
\alias{AnnotationIdentifier}
\alias{GOAllFrameIdentifier}
\alias{KEGGFrameIdentifier}
  
\title{Gene Identifier Class Constructors}
\description{
  Gene identifier classes and functions are used to indicate what the
  list of genes in a gene set represents (e.g., Entrez gene identifiers
  are tagged with \code{EntrezIdentifier()}, Bioconductor annotations with
  \code{AnnotationIdentifier()}).
}
\usage{
NullIdentifier(annotation, ...)
EnzymeIdentifier(annotation, ...)
ENSEMBLIdentifier(annotation, ...)
GenenameIdentifier(annotation,...)
RefseqIdentifier(annotation,...)
SymbolIdentifier(annotation,...)
UnigeneIdentifier(annotation,...)
EntrezIdentifier(annotation,...)
AnnotationIdentifier(annotation, ...)
}
\arguments{
  \item{annotation}{An optional character string identifying  the
    Bioconductor package from which the annotations are drawn, e.g.,
    \sQuote{hgu95av2}, \sQuote{org.Hs.eg.db}.}
  \item{\dots}{Additional arguments, usually none.}
}
\value{
  An object of the same class as the function name, initialized as
  appropriate for the identifier.
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{
  \code{\linkS4class{GeneIdentifierType}}-class for a description of the
  classes and methods using these objects.
}
\examples{
NullIdentifier()

data(sample.ExpressionSet)
gs1 <- GeneSet(sample.ExpressionSet[100:109],
               setName="sample1", setIdentifier="100")
geneIdType(gs1) # AnnotationIdentifier

geneIds <- featureNames(sample.ExpressionSet)[100:109]
gs2 <- GeneSet(geneIds=geneIds,
               setName="sample1", setIdentifier="101")
geneIdType(gs2) # NullIdentifier, since no info about genes provided

## Convert...
ai <- AnnotationIdentifier(annotation(sample.ExpressionSet))
geneIdType(gs2) <- ai
geneIdType(gs2)
## ...or provide more explicit construction
gs3 <- GeneSet(geneIds=geneIds, type=ai,
               setName="sample1", setIdentifier="102")
}
\keyword{manip}
