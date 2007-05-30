\name{GeneSet}
\docType{methods}
\alias{GeneSet}
\alias{GeneSet-methods}
\alias{GeneSet,missing-method}
\alias{GeneSet,character-method}
\alias{GeneSet,GeneIdentifierType-method}
\alias{GeneSet,ExpressionSet-method}
\title{Methods to construct GeneSet instances}
\description{
  Use \code{GeneSet} to construct gene sets from \code{ExpressionSet},
  character vector, or other objects.
}
\usage{
GeneSet(type, ...)
}
\arguments{
  \item{type}{An argument determining how the gene set will be created,
    as described in the Methods section.}
  \item{...}{Additional arguments for gene set construction. 
    Methods have required arguments, as outlined below; additional
    arguments correspond to slot names \code{\linkS4class{GeneSet}}.}
}
\section{Methods}{
  \describe{
    \item{\code{signature(type = "missing", setName =
        "character", setIdentifier = "character", ...)}}{Construct an
      empty gene set. The arguments \code{setName} and \code{setIdentifier} are
      required. The gene set has \code{setType=NullIdentifier},
      \code{collectionType=AdHocCollection}.}

    \item{\code{signature(type = "character", setName =
        "character", setIdentifier = "character", ...)}}{Construct an
      empty gene set. The arguments \code{setName} and \code{setIdentifier} are
      required. The gene set has \code{setType} created from the
      \code{GeneIdentifierType} named by \code{type}.}

    \item{\code{signature=(type = "GeneIdentifierType", setName =
        "character", setIdentifier = "character", ...)}}{Construct a
      empty gene set. The arguments \code{setName} and \code{setType} are
      required. The gene set has \code{setType} created from the
      \code{GeneIdentifierType} named by \code{type}.}

  \item{\code{signature=(type = "ExpressionSet", setName =
        "character", setIdentifier = "character", ...)}}{Construct a
      gene set from an
      \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. \code{setType} is set to
      \code{AnnotationIdentifier}; the annotation field and annotation package
      of the \code{ExpressionSet} are consulted to determine \code{organism},
      if possible.}
  }
}
\seealso{}
\examples{
## Empty gene set
GeneSet(setName="name", setIdentifier="100")
## Gene set from ExpressionSet
data(sample.ExpressionSet)
gs1 <- GeneSet(sample.ExpressionSet[100:109],
               setName="sample.GeneSet", setIdentifier="101")
## GeneSet from Broad XML; 'fl' could be a url
fl <- system.file("extdata", "Broad.xml", package="GSEABase")
gs2 <- getBroadSets(fl)[[1]] # actually, a list of two gene sets
## GeneSet from list of genes
genes <- genes(gs2) # any character vector would do
gs3 <- GeneSet(genes=genes, setName="another", setIdentifier="101")
## unspecified set type, so...
is(setType(gs3), "NullIdentifier") == TRUE
## update set type to match encoding of identifiers
setType(gs2)
setType(gs3) <- EntrezIdentifier()
}
\keyword{methods}
