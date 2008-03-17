\name{mapIdentifiers-methods}
\docType{methods}
\alias{mapIdentifiers}
\alias{mapIdentifiers-methods}
\alias{mapIdentifiers,ANY,ANY,ANY-method}
\alias{mapIdentifiers,GeneSet,GeneIdentifierType,missing-method}
\alias{mapIdentifiers,GeneSet,GeneIdentifierType,NullIdentifier-method}
\alias{mapIdentifiers,GeneSet,NullIdentifier,GeneIdentifierType-method}
\alias{mapIdentifiers,GeneSet,GeneIdentifierType,GeneIdentifierType-method}
\alias{mapIdentifiers,GeneSet,GeneIdentifierType,AnnDbBimap-method}
\alias{mapIdentifiers,GeneSet,GeneIdentifierType,environment-method}
\alias{mapIdentifiers,GeneSetCollection,GeneIdentifierType,missing-method}
\title{Methods for Function mapIdentifiers in Package `GSEABase'}
\description{
  These methods convert the genes identifiers of a gene set from one
  type to another, e.g., from \code{\link{EntrezIdentifier}} to
  \code{\link{AnnotationIdentifier}}. Methods can be called directly by
  the user; \code{\link{geneIdType<-}} provides similar
  functionality. \code{verbose=TRUE}  produces warning messages when
  maps between identifier types are not 1:1, or a map has to be
  constructed on the fly (this situation does not apply when using the
  DBI-based annotation packages).
}
\section{Methods}{
  \describe{

    \item{what = "ANY", to = "ANY", from = "ANY", verbose=FALSE}{
      This method warns of attempts to map \code{from} and \code{to} the
      same type, or generates an error if no suitable
      \code{mapIdentifiers} methods are available.}

    \item{what = "GeneSet", to = "GeneIdentifierType", from = "missing", verbose=FALSE}{
      This method will re-dispatch to a method with signature
      \code{signature(what=what, to=to, from=geneIdType(what))}, and is
      present so that a user can call \code{mapIdentifiers} without providing an
      explicit argument to \code{from}.}

    \item{what = "GeneSet", to = "GeneIdentifierType", from = "NullIdentifier", verbose=FALSE}{
      This maps a gene set from gene identifiers represented by the
      \code{NullIdentifier} type (i.e., no type associated with the genes) to
      gene identifiers represent by any class derived from
      \code{GeneIdentifierType}.}

    \item{what = "GeneSet", to = "NullIdentifier", from = "GeneIdentifierType", verbose=FALSE}{ 
      This maps a gene set from gene identifiers represented by any
      \code{GeneIdentifierType} type to one represented by the
      \code{NullIdentifier} (i.e., no type associated with the genes).}

    \item{what = "GeneSet", to = "GeneIdentifierType", from = "environment", verbose=FALSE}{
      Maps identifiers found in \code{what} to the type described by
      \code{to}, using the map (key-value pairs) found in \code{from}.}

    \item{what = "GeneSet", to = "GeneIdentifierType", from = "AnnDbBimap", verbose=FALSE}{
      Maps identifiers found in \code{what} to the type described by
      \code{to}, using the map (key-value pairs) found in \code{from}.}

    \item{what = "GeneSetCollection", to = "GeneIdentifierType", from = "missing", verbose = FALSE}{
      Map each gene set in \code{what} to gene identifier type
      \code{to}, using methods described above.}
}}
\keyword{methods}

