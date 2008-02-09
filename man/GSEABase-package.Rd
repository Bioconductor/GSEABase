\name{GSEABase-package}
\alias{GSEABase-package}
\alias{GSEABase}
\docType{package}
\title{
  Gene set enrichment data structures and methods
}
\description{
  This package provides classes and methods to support Gene Set
  Enrichment Analysis (GSEA).  The GeneSet class provides a
  common data structure for representing gene sets. The GeneColorSet
  class allows genes in a set to be associated with phenotypes. The
  GeneSetCollection class facilitates grouping together a list of related
  gene sets. The GeneIdentifierType class hierarchy reflects how genes
  are represented (e.g., Entrez versus symbol) in the gene
  set. mapIdentifiers provides a way to convert identifiers in a set from one
  type to another. The CollectionType class hierarchy reflects how the
  gene set was made, and can order genes into distinct sets or
  collections.
}
\author{
  Written by Martin Morgan, Seth Falcon, Robert Gentleman. 
  Maintainer: Biocore Team c/o BioC user list <bioconductor@stat.math.ethz.ch>
}
\seealso{
  \linkS4class{GeneSet},
  \linkS4class{GeneColorSet}
  \linkS4class{GeneSetCollection}
}
\examples{
example(GeneSet)
}
\keyword{package}
