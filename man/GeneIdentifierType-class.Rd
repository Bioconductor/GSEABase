\name{GeneIdentifierType-class}
\docType{class}
% GeneIdentifierType-class
\alias{GeneIdentifierType-class}
\alias{setType,GeneIdentifierType-method}
\alias{mapIdentifiers,NullIdentifier,GeneIdentifierType,GeneSet-method}
\alias{show,GeneIdentifierType-method}
\alias{validIdentifiers,GeneIdentifierType-method}
% AnnotationIdentifier-class
\alias{AnnotationIdentifier-class}
\alias{annotation,AnnotationIdentifier-method}
\alias{initialize,AnnotationIdentifier-method}
\alias{mapIdentifiers,AnnotationIdentifier,GeneIdentifierType,GeneSet-method}
\alias{mapIdentifiers,AnnotationIdentifier,character,GeneSet-method}
\alias{mapIdentifiers,AnnotationIdentifier,EntrezIdentifier,GeneSet-method}
\alias{show,AnnotationIdentifier-method}
% additional classes
\alias{NullIdentifier-class}
\alias{EntrezIdentifier-class}
\alias{PfamIdentifier-class}
\alias{SymbolIdentifier-class}

\title{Class "GeneIdentifierType"}

\description{
  This class provides a way to tag the meaning of gene
  symbols in a \code{GeneSet}. For instance, a \code{GeneSet} with gene
  names derived from a Bioconductor \code{annotation} package (e.g., via
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}) initially have a
  \code{GeneIdentifierType} of \code{AnnotationType}.
}

\section{Objects from the Class}{
  The following classes are available:
  \describe{
    \item{NullIdentifier}{No information about what genes represent.}
    \item{AnnotationIdentifier}{'Genes' are represented by the
      Affymetrix chip-specific probe identifier, as represented in Bicondcutor
      annotation packages.}
    \item{EntrezIdentifier}{Genes are represented by their Entrez
      identifier.}
    \item{PfamIdentifier}{Genes are represented by their Pfam
      identifier.}
    \item{SymbolIdentifier}{Genes are represented by their Symbol
      identifier.}
    \item{GeneIdentifierType}{A virtual Class: No objects may be created
      from it; all classes listed above are subclasses of
      \code{GeneIdentifierType}.}
  }
}

\section{Slots}{
  All \code{GeneIdentifierType} classes have the following slot:
  \describe{
    \item{type}{Object of class \code{"ScalarCharacter"}
      containing the character string representation of this
      \code{GeneIdentifierType}.}
  }

  The \code{AnnotationIdentifier} class has an additional slot
  \describe{
    \item{annotation}{Object of class \code{"ScalarCharacter"}
      containing the name of the annotation package from which the
      annotations (probe identifiers) are derived.}
  }
}
\section{Methods}{
  \code{GeneIdentifierType} classes are used in:
  \describe{
    \item{GeneSet}{\code{signature(type = "GeneIdentifierType")}:
      Create a new \code{\link{GeneSet}} using identifiers of
      \code{GeneIdentifierType}.}

    \item{GeneColorSet}{\code{signature(type = "GeneIdentifierType")}:
      Create a new \code{\link{GeneColorSet}} using identifiers of
      \code{GeneIdentifierType}.}

    \item{setType<-}{\code{signature(object = "GeneSet", value =
        "GeneIdentifierType")}: Change the \code{GeneIdentifierType} of
      \code{object} to \code{value}, attempting to convert symbols in
      the process.}

    \item{setType}{\code{signature(object = "GeneIdentifierType")}:
      return the type of this \code{object}.}

    \item{mapIdentifiers}{\code{\signature(from="NullIdentifier",
        to="GeneIdentifierType", whate="GeneSet")}: Used primarily
      internally (by \code{setType<-}), applied to a \code{geneSet} to map
      the set from genes represented by a \code{NullIdentifier} to genes
      represent by a class derived from \code{GeneIdentifierType}. Other
      \code{mapIdentifier} methods may change the representation of genes,
      but this method does not.}

  \item{show}{\code{signature(object = "GeneIdentifierType")}: display
      this object.}

    \item{validIdentifiers}{\code{signature(identifier =
        "GeneIdentifierType")}: validate that the identiifier are
      legitimate; subclasses are responsible for implementing this
      method.}
  }

  \code{AnnotationIdentifier} classes are also used in:
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
   The example below lists \code{GeneIdentifierType} classes defined in
   this package; See the help pages of these classes for specific information.
}
\examples{
getSubclasses(getClass("GeneIdentifierType"))

# create an AnnotationIdentifier, and ask it's type
setType(AnnotationIdentifier(annotation="hgu95av2"))

# Construct a GeneSet from an ExpressionSet, using the 'annotation'
# field of ExpressionSet to recognize the genes as AnnotationType
data(sample.ExpressionSet)
gs <- GeneSet(sample.ExpressionSet[100:109],
              setName="sample.GeneSet", setIdentifier="123")
setType(gs) # AnnotationIdentifier

## Read a Broad set from the sytem (or a url), and discover their
## GeneIdentifierType
fl <- system.file("extdata", "Broad.xml", package="GSEABase")
bsets <- getBroadSets(fl)
sapply(bsets, setType)

## try to combine gene sets with different set types
try(gs & sets[[1]])

\dontrun{
## Use the annotation package associated with the original
## ExpressionSet to map to EntrezIdentifier() ...
setType(gs) <- EntrezIdentifier()
...## and try again
gs & bsets[[1]]

## Another way to change annotation to Entrez (or other) ids
probeIds <- featureNames(sample.ExpressionSet)[100:109]
genes <- as.character(unlist(mget(probeIds, hgu95av2ENTREZID)))
GeneSet(EntrezIdentifier(),
        setName="sample.GenSet2", setIdentifier="101",
        genes=genes)
}

## Create a new identifier
setClass("FooIdentifier",
         contains="GeneIdentifierType",
         prototype=prototype(
           type=new("ScalarCharacter", "Foo")))
## Create a constructor (optional)
FooIdentifier <- function() new("FooIdentifier")
setType(FooIdentifier())

## tidy up
removeClass("FooIdentifier")
}
\keyword{classes}
