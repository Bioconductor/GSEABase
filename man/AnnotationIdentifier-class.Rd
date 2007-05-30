\name{AnnotationIdentifier-class}
\docType{class}
\alias{AnnotationIdentifier-class}
\alias{annotation,AnnotationIdentifier-method}
\alias{initialize,AnnotationIdentifier-method}
\alias{mapIdentifiers,AnnotationIdentifier,GeneIdentifierType,GeneSet-method}
\alias{mapIdentifiers,AnnotationIdentifier,character,GeneSet-method}
\alias{mapIdentifiers,AnnotationIdentifier,EntrezIdentifier,GeneSet-method}
\alias{show,AnnotationIdentifier-method}

\title{Class "AnnotationIdentifier"}
\description{A \code{\linkS4class{GeneIdentifierType}} class used to tag
  gene sets whose type of element is an Affymetrix annotation
  identifier, as represented by a \code{Bioconductor} annotation package.}
\section{Objects from the Class}{
  Objects can be created with \code{\link{AnnotationIdentifier}}; usually this
  is done during a call to \code{\link{GeneSet}} or \code{\link{GeneColorSet}}
}
\section{Slots}{
  \describe{
    \item{\code{annotation}:}{Object of class \code{"ScalarCharacter"}
      containing the string character representation of the chip from
      which the gene set was derived.}
    \item{\code{type}:}{Object of class \code{"ScalarCharacter"}
      containing a string identifier of this class.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{GeneIdentifierType}}, directly.
}
\section{Methods}{
  \describe{
    \item{annotation}{\code{signature(object = "AnnotationIdentifier")}:
      extract the name of the annotation package as a character string.}

    \item{mapIdentifiers}{\code{signature(from = "AnnotationIdentifier",
        to = "GeneIdentifierType", what = "GeneSet")}: Used primarily
      internally (by \code{setType<-}, applied to a \code{GeneSet}) to
      map between identifiers of type \code{AnnotationIdentifier} and other
      identifiers. This method consults the annotation package specified in
      \code{from}, looking for a mapping between \code{genes(what)} and
      \code{paste(annotation(what), toupper(setType(to)), sep="")}.}

    \item{mapIdentifiers}{\code{signature(from = "AnnotationIdentifier",
        to = "EntrezIdentifier", what = "GeneSet")}: Used primarily internally
      (by \code{setType<-}, applied to a \code{GeneSet}) to map between
      identifiers of type \code{AnnotationIdentifier} and
      \code{EntrezIdentifier}. This method consults the annotation
      package specified in \code{from}, looking for a mapping between
      \code{genes(what)} and \code{paste(annotation(what),
        "ENTREZID")}. }

    \item{mapIdentifiers}{\code{signature(from = "AnnotationIdentifier",
        to = "character", what = "GeneSet")}: Used primarily internally
      (by \code{setType<-}, applied to a \code{GeneSet}) to map between
      identifiers of type \code{AnnotationIdentifier} and other
      identifiers. This method differs from that with
      \code{to="GeneIdentifierType"} by allowing the user to specify an
      arbitrary character string for the mapping. The \code{setType} of
      the resulting set is either the result of a call to \code{to}, or
      the \code{\link{NullIdentifier}}.}

    \item{show}{\code{signature(object = "AnnotationIdentifier")}:
      Display \code{object}, including the annotation string if available.}
    \item{initialize}{\code{signature(.Object =
        "AnnotationIdentifier")}: Used internally during object creation.}
  }
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{
  \code{\linkS4class{GeneIdentifierType}}, \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}.
}
\examples{
AnnotationIdentifier()
example(GeneIdentifierType)
}
\keyword{classes}
