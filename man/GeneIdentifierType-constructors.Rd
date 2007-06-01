\name{GeneIdentifierType}
\alias{GeneIdentifierType}
\alias{NullIdentifier}
\alias{AnnotationIdentifier}
\alias{EntrezIdentifier}
\alias{SymbolIdentifier}
\alias{PfamIdentifier}
\title{Gene Identifier Class Constructors}
\description{
  Gene identifier classes and functions are used to indicate what the
  list of genes in a gene set represents (e.g., Entrez gene identifiers
  are tagged with \code{EntrezIdentifier()}, Bioconductor annotations with
  \code{AnnotationIdentifier()}).
}
\usage{
NullIdentifier(...)
AnnotationIdentifier(annotation,...)
EntrezIdentifier(...)
SymbolIdentifier(...)
PfamIdentifier(...)
}
\arguments{
  \item{annotation}{(required) A single character string identifying  the
    Bioconductor package from which the annotations are drawn.}
  \item{\dots}{Additional arguments, usually none.}
}
\value{
  An object of the same class as the function name, initialized as
  appropriate for the identifier.
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{
  \code{\linkS4class{GeneIdentifierType}} for a description of the
  classes and methods using these objects.
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
ai <- AnnotationIdentifier(annotation(sample.ExpressionSet))
setType(gs2) <- ai
setType(gs2)
## ...or provide more explicit construction
gs3 <- GeneSet(genes=genes, type=ai,
               setName="sample1", setIdentifier="102")
}
\keyword{manip}
