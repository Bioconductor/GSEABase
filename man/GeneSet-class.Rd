\name{GeneSet-class}
\docType{class}
\alias{GeneSet-class}
\alias{[,GeneSet,character,ANY-method}
\alias{[,GeneSet,numeric,ANY-method}
\alias{[,ExpressionSet,GeneSet,ANY-method}
\alias{[[,GeneSet,character-method}
\alias{[[,GeneSet,numeric-method}
\alias{Logic,character,GeneSet-method}
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
\alias{geneIds<-,GeneSet,character-method}
\alias{geneIds,GeneSet-method}
\alias{initialize,GeneSet-method}
\alias{intersect,GeneSet,GeneSet-method}
\alias{longDescription<-,GeneSet,character-method}
\alias{longDescription,GeneSet-method}
\alias{organism<-,GeneSet,character-method}
\alias{organism,GeneSet-method}
\alias{organism,GeneIdentifierType-method}
\alias{organism,GOAllFrameIdentifier-method}
\alias{organism,KEGGFrameIdentifier-method}
\alias{pubMedIds<-,GeneSet,character-method}
\alias{pubMedIds,GeneSet-method}
\alias{setdiff,GeneSet,GeneSet-method}
\alias{setIdentifier<-,GeneSet,character-method}
\alias{setIdentifier,GeneSet-method}
\alias{setName<-,GeneSet,character-method}
\alias{setName,GeneSet-method}
\alias{geneIdType<-,GeneSet,character-method}
\alias{geneIdType<-,GeneSet,GeneIdentifierType-method}
\alias{geneIdType,GeneSet-method}
\alias{setVersion<-,GeneSet,Versions-method}
\alias{setVersion,GeneSet-method}
\alias{show,GeneSet-method}
\alias{toGmt,GeneSet-method}
\alias{union,GeneSet,GeneSet-method}
\alias{urls<-,GeneSet,character-method}
\alias{urls,GeneSet-method}

\alias{collectionType<-}
\alias{collectionType}
\alias{contributor<-}
\alias{contributor}
\alias{creationDate<-}
\alias{creationDate}
\alias{geneIds<-}
\alias{geneIds}
\alias{geneIdType}
\alias{geneIdType<-}
\alias{longDescription<-}
\alias{longDescription}
\alias{organism<-}
\alias{organism}
\alias{setIdentifier<-}
\alias{setIdentifier}
\alias{setName<-}
\alias{setName}
\alias{setVersion<-}
\alias{setVersion}
\alias{urls<-}
\alias{urls}

\title{Class "GeneSet"}
\description{
  A \code{GeneSet} contains a set of gene identifiers. Each gene set has a
  \code{geneIdType}, indicating how the gene identifiers should be interpreted
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
}
\section{Slots}{
  \describe{
    \item{\code{setName}:}{Object of class \code{"ScalarCharacter"}
      containing a short name (single word is best) to identify the set.}
    \item{\code{setIdentifier}:}{Object of class
      \code{"ScalarCharacter"} containing a (unique) identifier for the
      set.}
    \item{\code{geneIdType}:}{Object of class \code{"GeneIdentifierType"}
      containing information about how the gene identifiers are encoded. See
      \code{\linkS4class{GeneIdentifierType}} and related classes.}
    \item{\code{geneIds}:}{Object of class \code{"character"} containing
      the gene symbols.}
    \item{\code{collectionType}:}{Object of class
      \code{"CollectionType"} containing information about how the geneIds
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
      containing the character string representation of the date on which
      the gene set was created.}
  }
}
\section{Methods}{
  Gene set construction:
  \describe{
    \item{GeneSet}{See \code{\link{GeneSet}} methods and
      \code{\link{getBroadSets}} for convenient construction.}
  }
  
  Slot access (e.g., \code{setName}) and retrieve
  (e.g., \code{setName<-}) :
  \describe{
    \item{collectionType<-}{\code{signature(object = "GeneSet", value = "CollectionType")}}
    \item{collectionType}{\code{signature(object = "GeneSet")}}
    \item{contributor<-}{\code{signature(object = "GeneSet", value = "character")}}
    \item{contributor}{\code{signature(object = "GeneSet")}}
    \item{creationDate<-}{\code{signature(object = "GeneSet", value = "character")}}
    \item{creationDate}{\code{signature(object = "GeneSet")}}
    \item{description<-}{\code{signature(object = "GeneSet", value = "character")}}
    \item{description}{\code{signature(object = "GeneSet")}}
    \item{geneIds<-}{\code{signature(object = "GeneSet", value = "character")}}
    \item{geneIds}{\code{signature(object = "GeneSet")}}
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
    \item{geneIdType<-}{
      \code{signature(object = "GeneSet", verbose=FALSE, value = "character")},
      \code{signature(object = "GeneSet", verbose=FALSE, value = "GeneIdentifierType")}:
      These method attempt to coerce geneIds from
      the current type to the type named by \code{value}. Successful
      coercion requires an appropriate method for \code{\link{mapIdentifiers}}.}
    \item{geneIdType}{\code{signature(object = "GeneSet")}}
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
      the same \code{geneIdType}.}
    \item{|}{
      \code{signature(e1 = "GeneSet", e2 = "character")},
      \code{signature(e1 = "character", e2 = "GeneSet")}:
      calculate the logical `or' (union) of a gene set and a character vector,
      i.e., add the geneIds named in the character vector to the gene set.}
    \item{intersect}{\code{signature(x = "GeneSet", y = "GeneSet")}:}
    \item{&}{\code{signature(e1 = "GeneSet", e2 = "GeneSet")}: calculate
      the logical `and' (intersection) of two gene sets.}
    \item{&}{
      \code{signature(e1 = "GeneSet", e2 = "character")},
      \code{signature(e1 = "character", e2 = "GeneSet")}:
      calculate the logical `and' (intersection) of a gene set and a
      character vector, creating a new gene set containing only those
      genes named in the character vector.}
    \item{setdiff}{
      \code{signature(x = "GeneSet", y = "GeneSet")},
      \code{signature(x = "GeneSet", y = "character")},
      \code{signature(x = "character", y = "GeneSet")}:
      calculate the logical set difference betwen two gene sets, or
      betwen a gene set and a character vector.}

    \item{[}{
      \code{signature(x = "GeneSet", i="character")}
      \code{signature(x = "GeneSet", i="numeric")}: subset the gene set by
      index (\code{i="numeric"}) or value (\code{i="character"}). Genes
      are re-ordered as required}
    \item{[}{
      \code{signature(x = "ExpressionSet", i = "GeneSet")}: subset the
      expression set, using genes in the gene set to select
      features. Genes in the gene set are coerced to appropriate annotation type
      if necessary (by consulting the \code{annotation} slot of the
      expression set, and using \code{geneIdType<-}).}
    \item{[[}{\code{signature(x = "GeneSet")}: select a single gene from
      the gene set.}
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
      methods. These convert genes from one \code{GeneIdentifierType} to another. See
      \code{\link{mapIdentifiers}} and specific methods in
      \code{\linkS4class{GeneIdentifierType}} for additional detail.}
    \item{incidence}{Summarize shared membership in genes across gene
      sets. See \code{\link{incidence-methods}}.}
    \item{toGmt}{Export to 'GMT' format file. See \code{\link{toGmt}}.}
    \item{show}{\code{signature(object = "GeneSet")}: display a short
      summary of the gene set.}
    \item{details}{\code{signature(object = "GeneSet")}: display
      additional information about the gene set. See \code{\link{details}}.}
    \item{initialize}{\code{signature(.Object = "GeneSet")}: Used
      internally during gene set construction.}
  }
}
\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{
  \code{\linkS4class{GeneColorSet}}
  \code{\linkS4class{CollectionType}}
  \code{\linkS4class{GeneIdentifierType}}
}
\examples{
## Empty gene set
GeneSet()
## Gene set from ExpressionSet
data(sample.ExpressionSet)
gs1 <- GeneSet(sample.ExpressionSet[100:109])
## GeneSet from Broad XML; 'fl' could be a url
fl <- system.file("extdata", "Broad.xml", package="GSEABase")
gs2 <- getBroadSets(fl)[[1]] # actually, a list of two gene sets
## GeneSet from list of geneIds
geneIds <- geneIds(gs2) # any character vector would do
gs3 <- GeneSet(geneIds=geneIds)
## unspecified set type, so...
is(geneIdType(gs3), "NullIdentifier") == TRUE
## update set type to match encoding of identifiers
geneIdType(gs2)
geneIdType(gs3) <- SymbolIdentifier()

## Convert between set types; this consults the 'annotation'
## information encoded in the 'AnnotationIdentifier' set type and the
## corresponding annotation package.
\dontrun{
gs4 <- gs1
geneIdType(gs4) <- EntrezIdentifier()
}

## logical (set) operations
gs5 <- GeneSet(sample.ExpressionSet[100:109], setName="subset1")
gs6 <- GeneSet(sample.ExpressionSet[105:114], setName="subset2")
## intersection: 5 'genes'; note the set name '(subset1 & subset2)'
gs5 & gs6
## union: 15 'genes'; note the set name
gs5 | gs6
## an identity
gs7 <- gs5 | gs6
gs8 <- setdiff(gs5, gs6) | (gs5 & gs6) | setdiff(gs6, gs5)
identical(geneIds(gs7), geneIds(gs8))
identical(gs7, gs8) == FALSE # gs7 and gs8 setNames differ

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
## list of geneIds --> vector of Broad GENESET XML
gs10 <- getBroadSets(fl) # two sets
entries <- sapply(gs10, function(x) toBroadXML(x))

## list mapIdentifiers available for GeneSet
showMethods("mapIdentifiers", classes="GeneSet", inherit=FALSE)
}
\keyword{classes}
