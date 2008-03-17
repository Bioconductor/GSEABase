\name{goSlim-methods}
\docType{methods}
\alias{goSlim}
\alias{goSlim-methods}
\alias{goSlim,GOCollection,GOCollection-method}
\alias{goSlim,ExpressionSet,GOCollection-method}
\title{Methods for Function goSlim in Package `GSEABase'}
\description{
  These methods summarize the gene ontology terms implied by the
  \code{idSrc} argument into the GO terms implied by the
  \code{slimCollection} argument. The summary takes identifiers in
  \code{idSrc} and determines all GO terms that
  apply to the identifiers. This full list of GO terms are then
  classified for membership in each term in the
  \code{slimCollection}.

  The resulting object is a data frame containing the terms of
  \code{slimCollection} as row labels, counts and frequencies of
  identifiers calssified to each term, and an abbreviated term
  description.

  An identifier in \code{idSrc} can expand to serveral GO terms, and the
  GO terms in \code{slimCollection} can imply an overlapping hierarchy
  of terms. Thus the resulting summary can easily contain more counts
  than there are identifiers in \code{idSrc}.
}
\usage{
goSlim(idSrc, slimCollection, ontology, ..., verbose=FALSE)
}
\arguments{
  \item{idSrc}{An argument determining the source of GO terms to be
    mapped to slim terms. The source might be a \code{GOCollection} of
    terms, or another object (e.g., ExpressionSet) for which the method
    can extract GO terms.}
  \item{slimCollection}{An argument containing the GO slim terms.}
  \item{ontology}{A character string nameing the ontology to be
    consulted when identifing slim term hierarchies. One of \sQuote{MF}
    (molecular function), \sQuote{BP} (biological process), \sQuote{CC}
    (cellular compartment).}
  \item{...}{Additional arguments passed to specific emthods.}
  \item{verbose}{Logical influencing whether messages (primarily missing
    GO terms arising during creation of the slim hierarchy) are reported.}
}
\section{Methods}{
  \describe{
    \item{idSrc="GOCollection", slimCollection="GOCollection",
      ontology="character", ..., verbose=FALSE}{Classify \code{idSrc} GO
      terms into \code{slimCollection} categories. The hierarchy of
      terms included for each term is from the ontology (MF, BP, or CC)
      specified by \code{ontology}. \code{verbose} informs about, e.g., GO
      terms that are not found.}
    \item{idSrc="ExpressionSet", slimCollection="GOCollection",
      ontology="character", ...,  verbose=FALSE}{Determine the (unique) GO terms
      implied by feature names in \code{idSrc} (using the annotation map
      identified in \code{annotation(idSrc)}).}
}}
\examples{

myIds <- c("GO:0016564", "GO:0003677", "GO:0004345", "GO:0008265",
           "GO:0003841", "GO:0030151", "GO:0006355", "GO:0009664",
           "GO:0006412", "GO:0015979", "GO:0006457", "GO:0005618",
           "GO:0005622", "GO:0005840", "GO:0015935", "GO:0000311",
           "GO:0009282")
myCollection <- GOCollection(myIds)
fl <- system.file("extdata", "goslim_plant.obo", package="GSEABase")
slim <- getOBOCollection(fl)
goSlim(myCollection, slim, "MF")
data(sample.ExpressionSet)
goSlim(sample.ExpressionSet, slim, "MF", evidenceCode="TAS")
}
\keyword{methods}

