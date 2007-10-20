\name{CollectionType-class}
\docType{class}
% CollectionType-class
\alias{CollectionType-class}
\alias{collectionType,CollectionType-method}
\alias{show,CollectionType-method}
\alias{CollectionIdType-class}
\alias{[,CollectionIdType,missing,missing-method}
\alias{show,CollectionIdType-method}
\alias{ids}
\alias{ids,CollectionIdType-method}
% NullCollection-class
\alias{NullCollection-class}
% ExpressionSetCollection-class
\alias{ExpressionSetCollection-class}
% BroadCollection-class
\alias{BroadCollection-class}
\alias{bcCategory,BroadCollection-method}
\alias{bcSubCategory,BroadCollection-method}
\alias{show,BroadCollection-method}
\alias{bcCategory}
\alias{bcSubCategory}
% GOCollection-class
\alias{GOCollection-class}
\alias{evidenceCode,GOCollection-method}
\alias{goIds,GOCollection-method}
\alias{[,GOCollection,missing,missing-method}
\alias{show,GOCollection-method}
\alias{evidenceCode}
\alias{goIds}
% KEGGCollection-class
\alias{KEGGCollection-class}
% additional classes
\alias{ChrCollection-class}
\alias{ChrlocCollection-class}
\alias{OMIMCollection-class}
\alias{PMIDCollection-class} 
\alias{PfamCollection-class}
\alias{PrositeCollection-class}

\title{Class "CollectionType"}

\description{
  These classes provides a way to tag the origin of a
  \code{GeneSet}. Collection types can be used in manipulating
  (e.g., selecting) sets, and can contain information specific to
  particular sets (e.g., \code{category} and \code{subcategory}
  classifications of \code{BroadCollection}.)
}

\section{Objects from the Class}{
  The following classes can tag gene sets; \code{GO, KEGG, Chr, Chrloc,
    OMIM, and PMID} collections can be derived from chip or organism
  \sQuote{annotation} packages.

  \describe{
    \item{NullCollection}{No formal collection information available.}
    \item{BroadCollection}{Derived from, or destined to be, Broad
      XML. Usually created and written \code{\link{getBroadSets}},
      \code{\link{toBroadXML}}.}
    \item{ExpressionSetCollection}{Derived from
      \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. Usually
      created during a call to \code{\link{GeneSet}} or
      \code{\link{GeneColorSet}}.}
    \item{GOCollection}{Collection derived using Gene Ontology (GO)
      terms.}
    \item{OBOCollection}{Collection derived from \code{GOCollection},
      specifically from files described by the OBO file format. See
      \code{\linkS4class{OBOCollection}}}
    \item{KEGGCollection}{Collection derived using KEGG terms.}
    \item{ChrCollection}{Collection derived using chromsome locations}
    \item{ChrlocCollection}{Collection derived using chromosome starting
      posistions}
    \item{MapCollection}{Collection derived from cytogenic bands.}
    \item{OMIMCollection}{Collection derived from identifiers in the
      Online Inheritance in Man.}
    \item{PMIDCollection}{Collection derived from PMID identifiers.}
    \item{PfamCollection}{Collection derived from Pfam identifiers.}
    \item{PrositeCollection}{Collection derived from Prosite identifiers.}
  }

  Objects are instantiated with calls to \code{\link{CollectionType}}
  constructors, with slot names as possible arguments.
}

\section{Slots}{
  \code{CollectionType} classes (\code{Null, ExpressionSet}) have the slot:
  \describe{
    \item{\code{type}:}{Object of class \code{"ScalarCharacter"}
      containing the character string representation of this
      \code{CollectionType}.}
  }

  \code{CollectionIdType} classes (\code{KEGG, OMIM, PMID, Chr, Chrloc,
    Map, GO}) extend the \code{CollectionType} and have the additional
  slot:
  \describe{
    \item{\code{ids}:}{Object of class \code{"character"} containing a
      vector of character string representations of corresponding
      identifiers, e.g., \sQuote{KEGG} or \sQuote{GO} terms.}
  }

  \code{GOCollection} extends \code{CollectionIdType} and has the
  additional slot:
  \describe{
    \item{\code{evidenceCode}:}{Object of class \code{"character"},
      containing GO evidence codes used to construct the gene set.}
  }

  \code{OBOCollection} extends \code{GOCollection}; see
  \code{\linkS4class{OBOCollection}}.

  \code{BroadCollection} has slots:
  \describe{
     \item{\code{category}:}{Object of class \code{"ScalarCharacter"}
      containing terms from the Broad list of categories, or \code{NA}}
    \item{\code{subCategory}:}{Object of class \code{"ScalarCharacter"}
      containing Broad sub-categories, or \code{NA}}
  }
}

\section{Methods}{
  \code{CollectionType} classes have methods:
  \describe{
    \item{collectionType<-}{\code{signature(object = "GeneSet", value =
        "CollectionType")}: Replace the \code{CollectionType} }
    \item{collectionType}{\code{signature(object = "CollectionType")}:
      Retrieve the collection type.}
    \item{show}{\code{signature(object = "CollectionType")}: display the
      collection type.}
  }

  \code{CollectionIdType} classes inherit \code{CollectionType} methods,
  and have in addition:
  \describe{
    \item{ids}{\code{signature(object="CollectionIdType")}: Retrieve the
      identifiers of the collection type.}
    \item{[}{\code{signature(object="CollectionIdType", i="missing",
        j="missing", ..., ids=ids(object))}: return a subset of
      \code{object} containing only ids in \code{ids}}
  }

  \code{GOCollection} inherits \code{CollectionIdType} methods, and has
  in addition:
  \describe{
    \item{evidenceCode}{Retrieve the evidence codes of the GO
      collection. Evidence codes are summarize in, e.g.,
      \code{\link[GO]{GOENTREZID}}.}
    \item{[}{\code{signature(object="CollectionIdType", i="missing",
        j="missing", ..., evidenceCode=evidenceCode(object))}: return a subset of
      \code{object} containing only evidence codes in
      \code{evidenceCode}. This method passes arguments ... to
      \code{[,CollectionIdType} methods.}
  }

  \code{BroadCollection} has methods:
  \describe{
    \item{bcCategory}{Retrieve the \code{category} of the Broad
      collection.}
    \item{bcSubCategory}{Retrieve the sub-category of the Broad
      collection.}
  }
}

\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{
  \code{\link{CollectionType}} consturctors.
}
\examples{
getSubclasses(getClass("CollectionType"))

## Create a CollectionType and ask for its type
collectionType(ExpressionSetCollection())

\dontrun{
## Read two GeneSets from a Broad XML file into a list, verify that
## they are both BroadCollection's. Category / subcategory information
## is unique to Broad collections.
fl <- system.file("extdata", "Broad.xml", package="GSEABase")
sets <- getBroadSets(fl)
sapply(sets, collectionType)
}

## ExpressionSets are tagged with ExpressionSetCollection; there is no
## 'category' information.
data(sample.ExpressionSet)
gs <- GeneSet(sample.ExpressionSet[100:109],
              setName="sample.GeneSet", setIdentifier="123")
collectionType(gs)

## GOCollections are created by reference to GO terms and evidenceCodes
GOCollection("GO:0005488")
## requires library(GO); EntrezIdentifers automatically created
\dontrun{
GeneSet(GOCollection(c("GO:0005488", "GO:0019825"),
                     evidenceCode="IDA"))
}
}
\keyword{classes}
