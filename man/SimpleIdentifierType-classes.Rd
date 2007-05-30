\name{Simple Identifiers}
\docType{class}
\alias{<name>Identifier-class}
\alias{NullIdentifier-class}
\alias{EntrezIdentifier-class}
\alias{PfamIdentifier-class}
\alias{SymbolIdentifier-class}

\title{Simple "GeneIdentifierType" Classes}
\description{

  These \code{\linkS4class{GeneIdentifierType}} classes are used to tag
  gene sets to indicate what gene elements represent.Available simple
  classes incldue:
  \describe{
    \item{NullIdentifier}{No information about what genes represent.}
    \item{EntrezIdentifier}{Genes are represented by their Entrez
      identifier.}
    \item{PfamIdentifier}{Genes are represented by their Pfam
      identifier.}
    \item{SymbolIdentifier}{Genes are represented by their Symbol
      identifier.}
  }

  Additional classes extending \code{GeneIdentifierType} are available,
  especialy \code{\linkS4class{AnnotationIdentifier}} for representing
  "genes" (i., probe sets) derived from \code{ExpressionSet}.
}
\section{Objects from the Class}{
  Objects can be created with \code{\link{NullIdentifier}} and similar
  functions names after the gene identier; usually this
  is done during a call to \code{\link{GeneSet}} or
  \code{\link{GeneColorSet}}.
}

\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"ScalarCharacter"}
      containing a string identifier of this class.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{GeneIdentifierType}}, directly.
}
\section{Methods}{
  Methods for \code{setType}, \code{setType<-}, and
  \code{mapIdentifiers} are defined on these types.
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{
  \code{\linkS4class{GeneIdentifierType}},
  \code{\link{setType}},
  \code{\link{setType<-}},
  \code{\link{mapIdentifiers}}.
}
\examples{
example(GeneIdentifierType)
}
\keyword{classes}
