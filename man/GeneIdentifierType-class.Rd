\name{GeneIdentifierType-class}
\docType{class}
% GeneIdentifierType-class
\alias{GeneIdentifierType-class}
\alias{geneIdType,GeneIdentifierType-method}
\alias{show,GeneIdentifierType-method}
% AnnotationIdentifier-class
\alias{AnnotationIdentifier-class}
\alias{annotation,AnnotationIdentifier-method}
\alias{initialize,AnnotationIdentifier-method}
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
  \code{GeneIdentifierType} of \code{AnnotationIdentifier}.
}

\section{Objects from the Class}{
  The following classes are available:
  \describe{
    \item{NullIdentifier}{No formal information about what gene identifiers represent.}
    \item{AnnotationIdentifier}{Gene identifiers are Affymetrix
      chip-specific probe identifier, as represented in Bioconductor
      annotation packages.}
    \item{EntrezIdentifier}{Gene identifiers \sQuote{Entrez} identifiers.}
    \item{PfamIdentifier}{Gene identifiers are \sQuote{Pfam} identifiers.}
    \item{SymbolIdentifier}{Genes identifiers are \sQuote{Symbol} identifiers.}
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

    \item{geneIdType<-}{\code{signature(object = "GeneSet", value =
        "GeneIdentifierType")}: Changes the \code{GeneIdentifierType} of
      \code{object} to \code{value}, attempting to convert symbols in
      the process. This method calls \code{mapIdentifiers(what=object,
        to=value, from=geneIdType(what))}.}

    \item{geneIdType}{\code{signature(object = "GeneIdentifierType")}:
      return a character string representation of the type of this
      \code{object}.}

    \item{mapIdentifiers}{
      \code{signature(what="GeneSet", to="GeneIdentifierType" from="missing")},
      \code{signature(what="GeneSet", to="GeneIdentifierType" from="NullIdentifier")}
      : See \code{\link{mapIdentifiers}}.}

    \item{show}{\code{signature(object = "GeneIdentifierType")}: display
      this object.}
  }

  \code{AnnotationIdentifier} has the
  following method:
  \describe{
    \item{geneIdType}{\code{signature(object = "AnnotationIdentifier")}:
      Retrieve the \code{geneIdType} (e.g., annotation string for
      \code{AnnotationIdentifier}; ad hoc string for \code{NullIdentifier}).}
  }

  \code{AnnotationIdentifier} classes are also used in:
  \describe{
    \item{annotation}{\code{signature(object = "AnnotationIdentifier")}:
      extract the name of the annotation package as a character string.}

    \item{mapIdentifiers}{
      \code{signature(what = "GeneSet", to = "GeneIdentifierType", from = "AnnotationIdentifier")},
      \code{signature(what = "GeneSet", to="AnnotationIdentifier", from = "GeneIdentifierType")}:
      see \code{\link{mapIdentifiers}}.}

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
geneIdType(AnnotationIdentifier(annotation="hgu95av2"))

# Construct a GeneSet from an ExpressionSet, using the 'annotation'
# field of ExpressionSet to recognize the genes as AnnotationType
data(sample.ExpressionSet)
gs <- GeneSet(sample.ExpressionSet[100:109],
              setName="sample.GeneSet", setIdentifier="123")
geneIdType(gs) # AnnotationIdentifier

## Read a Broad set from the system (or a url), and discover their
## GeneIdentifierType
fl <- system.file("extdata", "Broad.xml", package="GSEABase")
bsets <- getBroadSets(fl)
sapply(bsets, geneIdType)

## try to combine gene sets with different set types
try(gs & sets[[1]])

\dontrun{
## Use the annotation package associated with the original
## ExpressionSet to map to EntrezIdentifier() ...
geneIdType(gs) <- EntrezIdentifier()
...## and try again
gs & bsets[[1]]

## Another way to change annotation to Entrez (or other) ids
probeIds <- featureNames(sample.ExpressionSet)[100:109]
geneIds <- as.character(unlist(mget(probeIds, hgu95av2ENTREZID)))
GeneSet(EntrezIdentifier(),
        setName="sample.GenSet2", setIdentifier="101",
        geneIds=geneIds)
}

## Create a new identifier
setClass("FooIdentifier",
         contains="GeneIdentifierType",
         prototype=prototype(
           type=new("ScalarCharacter", "Foo")))
## Create a constructor (optional)
FooIdentifier <- function() new("FooIdentifier")
geneIdType(FooIdentifier())

## tidy up
removeClass("FooIdentifier")
}
\keyword{classes}
