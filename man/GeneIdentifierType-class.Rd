\name{GeneIdentifierType-class}
\docType{class}
% GeneIdentifierType-class
\alias{GeneIdentifierType-class}
\alias{geneIdType,GeneIdentifierType-method}
\alias{show,GeneIdentifierType-method}
% AnnotationIdentifier-class
\alias{AnnotationIdentifier-class}
\alias{annotation,GeneIdentifierType-method}
\alias{annotation<-,GeneIdentifierType,character-method}
\alias{initialize,GeneIdentifierType-method}
\alias{show,AnnotationIdentifier-method}
% additional classes
\alias{NullIdentifier-class}
\alias{EnzymeIdentifier-class}
\alias{GenenameIdentifier-class}
\alias{RefseqIdentifier-class}
\alias{SymbolIdentifier-class}
\alias{UnigeneIdentifier-class}
\alias{EntrezIdentifier-class}

\title{Class "GeneIdentifierType"}

\description{
  This class provides a way to tag the meaning of gene
  symbols in a \code{GeneSet}. For instance, a \code{GeneSet} with gene
  names derived from a Bioconductor \code{annotation} package (e.g., via
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}) initially have a
  \code{GeneIdentifierType} of \code{AnnotationIdentifier}.
}

\section{Objects from the Class}{
  The following classes are available, and derive from tables in 
  \sQuote{annotation} packages
  \describe{
    \item{NullIdentifier}{No formal information about what gene identifiers represent.}
    \item{AnnotationIdentifier}{Gene identifiers are Affymetrix
      chip-specific probe identifier, as represented in Bioconductor
      annotation packages.}
    \item{EntrezIdentifier}{\sQuote{Entrez} identifiers.}
    \item{EnzymeIdentifier}{\sQuote{EC} identifiers.}
    \item{GenenameIdentifier}{Curated and ad hoc descriptive gene names.}
    \item{RefseqIdentifier}{\sQuote{Prosite} identifiers.}
    \item{SymbolIdentifier}{\sQuote{Symbol} identifiers.}
    \item{UnigeneIdentifier}{\sQuote{Unigene} identifiers.}
    \item{GeneIdentifierType}{A virtual Class: No objects may be created
      from it; all classes listed above are subclasses of
      \code{GeneIdentifierType}.}
  }
}

\section{Slots}{
  All \code{GeneIdentifierType} classes have the following slots:
  \describe{
    \item{type}{Object of class \code{"ScalarCharacter"}
      containing the character string representation of this
      \code{GeneIdentifierType}.}
    \item{annotation}{Object of class \code{"ScalarCharacter"}
      containing the name of the annotation package from which the
      identifiers (probe identifiers) are derived.}
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

    \item{annotation}{\code{signature(object = "GeneIdentifierType")}:
      extract the name of the annotation package as a character string.}

    \item{annotation<-}{\code{signature(object = "GeneIdentifierType", value = "character")}:
      assign the name of the annotation package as a character string.}

    \item{geneIdType}{\code{signature(object = "GeneIdentifierType")}:
      return a character string representation of the type of this
      \code{object}.}

    \item{geneIdType<-}{\code{signature(object = "GeneSet", verbose=FALSE, value =
        "GeneIdentifierType")}: Changes the \code{GeneIdentifierType} of
      \code{object} to \code{value}, attempting to convert symbols in
      the process. This method calls \code{mapIdentifiers(what=object,
        to=value, from=geneIdType(what), verbose=verbose)}.}

    \item{mapIdentifiers}{See \code{\link{mapIdentifiers}}.}

    \item{show}{\code{signature(object = "GeneIdentifierType")}: display
      this object.}
  }
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{
   The example below lists \code{GeneIdentifierType} classes defined in
   this package; See the help pages of these classes for specific information.
}
\examples{
names(getClass("GeneIdentifierType")@subclasses)

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
## ...and try again
gs & bsets[[1]]

## Another way to change annotation to Entrez (or other) ids
probeIds <- featureNames(sample.ExpressionSet)[100:109]
geneIds <- getEG(probeIds, "hgu95av2")
GeneSet(EntrezIdentifier(),
        setName="sample.GeneSet2", setIdentifier="101",
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
