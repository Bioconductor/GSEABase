\name{GeneSet-class}
\docType{class}
\alias{GeneSet-class}
\alias{[,GeneSet,character-method}
\alias{[,GeneSet,numeric-method}
\alias{[[,GeneSet,character-method}
\alias{[[,GeneSet,numeric-method}
\alias{|,GeneSet,GeneSet-method}
\alias{|,GeneSet,character-method}
\alias{$,GeneSet-method}
\alias{&,GeneSet,GeneSet-method}
\alias{&,GeneSet,character-method}
\alias{collectionType<-,GeneSet,CollectionType-method}
\alias{collectionType,GeneSet-method}
\alias{contributor<-,GeneSet,character-method}
\alias{contributor,GeneSet-method}
\alias{creationDate<-,GeneSet,character-method}
\alias{creationDate,GeneSet-method}
\alias{description<-,GeneSet,character-method}
\alias{description,GeneSet-method}
\alias{genes<-,GeneSet,character-method}
\alias{genes,GeneSet-method}
\alias{initialize,GeneSet-method}
\alias{intersect,GeneSet,GeneSet-method}
\alias{longDescription<-,GeneSet,character-method}
\alias{longDescription,GeneSet-method}
\alias{organism<-,GeneSet,character-method}
\alias{organism,GeneSet-method}
\alias{pubMedIds<-,GeneSet,character-method}
\alias{pubMedIds,GeneSet-method}
\alias{setdiff,GeneSet,GeneSet-method}
\alias{setIdentifier<-,GeneSet,character-method}
\alias{setIdentifier,GeneSet-method}
\alias{setName<-,GeneSet,character-method}
\alias{setName,GeneSet-method}
\alias{setType<-,GeneSet,character-method}
\alias{setType<-,GeneSet,GeneIdentifierType-method}
\alias{setType,GeneSet-method}
\alias{setVersion<-,GeneSet,Versions-method}
\alias{setVersion,GeneSet-method}
\alias{show,GeneSet-method}
\alias{union,GeneSet,GeneSet-method}
\alias{urls<-,GeneSet,character-method}
\alias{urls,GeneSet-method}

\title{Class "GeneSet"}
\description{
  A \code{GeneSet} contains a set of genes names. Each gene set has a
  \code{setType}, indicating how the gene names should be interpreted
  (e.g., as Entrez identifiers), and a \code{collectionType}, indicating
  the origin of the gene set (perhaps including additional information
  about the set, as in the \code{\linkS4class{BroadCollection}} type).

  Conversion between identifiers, subsetting, and logical (set)
  operations can be performed. Relationships between genes and phenotype
  in a \code{GeneSet} can be summarized using \code{coloring} to create
  a \code{GeneColorSet}. A \code{GeneSet} can be exported to XML with
  \code{toBroadXML}. 
}
\section{Objects from the Class}{
  Construct a \code{GeneSet} with a \code{\link{GeneSet}} method (e.g.,
  from a character vector of gene names, or an
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}), or from gene sets stored as XML
  (locally or on the internet; see \code{getBroadSets})

  A \code{GeneSet} requires a \code{setName} (best just to be a single
  word) and a unique \code{setIdentifier}. Aditional informational
  elements can be added (see Slots, below, for an indication of this
  information, and the Methods for how to set or retrieve it.).
}
\section{Slots}{
	 \describe{
    \item{\code{setName}:}{Object of class \code{"ScalarCharacter"}
      containing a short name (single word is best) to identify the set.}
    \item{\code{setIdentifier}:}{Object of class
      \code{"ScalarCharacter"} containing a (unique) identifier for the
      set.}
    \item{\code{type}:}{Object of class \code{"GeneIdentifierType"}
      containing information about how the genes are encoded. See
      \code{\linkS4class{GeneIdentifierType}} and related classes.}
    \item{\code{genes}:}{Object of class \code{"character"} containing
      the gene symbols.}
    \item{\code{collectionType}:}{Object of class
      \code{"CollectionType"} containing information about how the genes
      were collected, including perhaps additional information unique to
      the collection methodology. See \code{\linkS4class{CollectionType}}
      and related classes.}
    \item{\code{shortDescription}:}{Object of class
      \code{"ScalarCharacter"} representing short description (1 line) of the gene set.}
    \item{\code{longDescription}:}{Object of class
      \code{"ScalarCharacter"} providing a longer description (e.g.,
      like an abstract) of the gene set.}
    \item{\code{organism}:}{Object of class \code{"ScalarCharacter"}
      represents the organism the gene set is derived from.}
    \item{\code{pubMedIds}:}{Object of class \code{"character"}
      containing PubMed ids related to the gene set.}
    \item{\code{urls}:}{Object of class \code{"character"} containing
      urls used to construct or manipulate the gene set.}
    \item{\code{contributor}:}{Object of class \code{"character"}
      identifying who created the gene set.}
    \item{\code{version}:}{Object of class \code{"Versions"} a version
      number, manually curated (i.e., by the \code{contributor}) to
      provide a consistent way of tracking a gene set.}
    \item{\code{creationDate}:}{Object of class \code{"character"}
      containg the character string representation of the date on which
      the gene set was created.}
  }
}
\section{Methods}{
  Gene set construction:
  \describe{
    \item{GeneSet}{See \code{\link{GeneSet}} methods and
      \code{\link{getBroadSets}} for convenient construction.}
  }
  
  Slot access (e.g., \code{collectionType}) and retrieve
  (e.g., \code{collectionType<-}) :
  \describe{
    \item{collectionType<-}{\code{signature(object = "GeneSet", value = "CollectionType")}}
    \item{collectionType}{\code{signature(object = "GeneSet")}}
    \item{contributor<-}{\code{signature(object = "GeneSet", value = "character")}}
    \item{contributor}{\code{signature(object = "GeneSet")}}
    \item{creationDate<-}{\code{signature(object = "GeneSet", value = "character")}}
    \item{creationDate}{\code{signature(object = "GeneSet")}}
    \item{description<-}{\code{signature(object = "GeneSet", value = "character")}}
    \item{description}{\code{signature(object = "GeneSet")}}
    \item{genes<-}{\code{signature(object = "GeneSet", value = "character")}}
    \item{genes}{\code{signature(object = "GeneSet")}}
    \item{longDescription<-}{\code{signature(object = "GeneSet", value = "character")}}
    \item{longDescription}{\code{signature(object = "GeneSet")}}
    \item{organism<-}{\code{signature(object = "GeneSet", value = "character")}}
    \item{organism}{\code{signature(object = "GeneSet")}}
    \item{pubMedIds<-}{\code{signature(object = "GeneSet", value = "character")}}
    \item{pubMedIds}{\code{signature(object = "GeneSet")}}
    \item{setdiff}{\code{signature(x = "GeneSet", y = "GeneSet")}}
    \item{setIdentifier<-}{\code{signature(object = "GeneSet", value = "character")}}
    \item{setIdentifier}{\code{signature(object = "GeneSet")}}
    \item{setName<-}{\code{signature(object = "GeneSet", value = "character")}}
    \item{setName}{\code{signature(object = "GeneSet")}}
    \item{setType<-}{
      \code{signature(object = "GeneSet", value = "character")},
      \code{signature(object = "GeneSet", value = "GeneIdentifierType")}:
      These method attempt to coerce genes from
      the current type to the type named by \code{value}. Successful
      coersion requires an appropariate method for \code{\link{mapIdentifiers}}.}
    \item{setType}{\code{signature(object = "GeneSet")}}
    \item{setVersion<-}{\code{signature(object = "GeneSet", value = "Versions")}}
    \item{setVersion}{\code{signature(object = "GeneSet")}}
    \item{urls<-}{\code{signature(object = "GeneSet", value = "character")}}
    \item{urls}{\code{signature(object = "GeneSet")}}
  }

  Logical and subsetting operations:
  \describe{
    \item{union}{\code{signature(x = "GeneSet", y = "GeneSet")}: ... }
    \item{|}{\code{signature(e1 = "GeneSet", e2 = "GeneSet")}: calculate
      the logical `or' (union) of two gene sets. The sets must contain elements of
      the same \code{setType}.}
    \item{|}{\code{signature(e1 = "GeneSet", e2 = "character")}:
      calculate the logcial `or' (union) of a gene set and a character vector,
      i.e., determine the genes named in the character vector also
      present in the gene set.}
    \item{intersect}{\code{signature(x = "GeneSet", y = "GeneSet")}:{}}
    \item{&}{\code{signature(e1 = "GeneSet", e2 = "GeneSet")}: calculate
    the logical `and' (intersection) of two gene sets.}
    \item{&}{\code{signature(e1 = "GeneSet", e2 = "character")}:
      calculate the logcial `and' (intersection) of a gene set and a
      character vector, creating a new gene set with the character
      vector coerced to be of type \code{setType}.}

    \item{[}{
      \code{signature(x = "GeneSet", i="character")}
      \code{signature(x = "GeneSet", i="numeric")}: subset the gene set by
      index (\code{i="numeric"}) or value (\code{i="character"}). Genes
      are re-ordered as required}
    \item{[[}{\code{signature(x = "GeneSet")}: select a single gene from
      the gene set; not particularly useful.}
    \item{\$}{\code{signature(x = "GeneSet")}: select a single gene from
      the gene set, allowing partial matching.}
  }

  Useful additional methods include:
  \describe{
    \item{GeneColorSet}{\code{signature(type = "GeneSet")}: create a
      'color' gene set from a \code{GeneSet}, containing information
      about phenotype. This method has a required argument
      \code{phenotype}, a character string describing the phenotype for
      which color is available. See \code{\linkS4class{GeneColorSet}}.}
    \item{mapIdentifiers}{Use the code in the examples to list available
      methods. These are used primarily internally (by \code{setType<-})
      to convert gene from one \code{GeneIdentifierType} to another. See
      \code{\link{mapIdentifiers}} for additional detail.}
    \item{show}{\code{signature(object = "GeneSet")}: display the gene set.}
    \item{initialize}{\code{signature(.Object = "GeneSet")}: Used
      internally during gene set construction.}
  }
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{
  \code{\linkS4class{CollectionType}}
  \code{\linkS4class{GeneIdentifierType}}
}
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

## Convert between set types; this consults the 'annotation'
## information encoded in the 'AnnotationIdentifier' set type and the
## correspondining annotation package.
\dontrun{
gs4 <- gs1
setType(gs4) <- EntrezIdentifier()
}

## logical (set) operations
gs5 <- GeneSet(sample.ExpressionSet[100:109],
               setName="subset1", setIdentifier="105")
gs6 <- GeneSet(sample.ExpressionSet[105:114],
               setName="subset2", setIdentifier="106")
## intersection: 5 'genes'; note the set name '(subset1 & subset2)'
gs5 & gs6
## union: 15 'genes'; note the set name
gs5 | gs6
## an identity
gs7 <- gs5 | gs6
gs8 <- setdiff(gs5, gs6) | (gs5 & gs6) | setdiff(gs6, gs5)
identical(genes(gs7), genes(gs8))
identical(gs7, gs8) == FALSE # gs7 and gs8 setNmes differ

## output
tmp <- tempfile()
toBroadXML(gs2, tmp)
noquote(readLines(tmp))
## must be BroadCollection() collectionType 
try(toBroadXML(gs1))
gs9 <- gs1
collectionType(gs9) <- BroadCollection()
toBroadXML(gs9, tmp)
unlink(tmp)
toBroadXML(gs9) # no connection --> character vector
## list of genes --> vector of Broad GENESET XML
gs10 <- getBroadSets(fl) # two sets
entries <- sapply(gs10, function(x) toBroadXML(x)[[2]])

## list mapIdentifiers available for GeneSet
showMethods("mapIdentifiers", classes="GeneSet", inherit=FALSE)
}
\keyword{classes}
