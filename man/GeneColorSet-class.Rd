\name{GeneColorSet-class}
\docType{class}
\alias{GeneColorSet-class}
\alias{|,GeneColorSet,GeneColorSet-method}
\alias{|,GeneColorSet,character-method}
\alias{&,GeneColorSet,GeneColorSet-method}
\alias{&,GeneColorSet,character-method}
\alias{coloring<-,GeneColorSet,data.frame-method}
\alias{coloring,GeneColorSet-method}
\alias{geneColor<-,GeneColorSet,factor-method}
\alias{geneColor,GeneColorSet-method}
\alias{initialize,GeneColorSet-method}
\alias{intersect,GeneColorSet,GeneColorSet-method}
\alias{phenotypeColor<-,GeneColorSet,factor-method}
\alias{phenotypeColor,GeneColorSet-method}
\alias{phenotype<-,GeneColorSet,character-method}
\alias{phenotype,GeneColorSet-method}
\alias{setdiff,GeneColorSet,GeneColorSet-method}
\alias{show,GeneColorSet-method}
\alias{union,GeneColorSet,GeneColorSet-method}
\alias{[,GeneColorSet,numeric-method}
\alias{[,GeneColorSet,character-method}
\alias{[[,GeneColorSet,numeric-method}
\alias{[[,GeneColorSet,character-method}
\alias{$,GeneColorSet-method}


\title{Class "GeneColorSet"}
\description{
  A \code{GeneColorSet} extends \code{\linkS4class{GeneSet}} to allow
  genes to be 'colored'. Coloring means that for a particular phenotype,
  each gene has a color (e.g., expression levels "up", "down", or "unchanged")
  and a phenotypic consequence (e.g., the phenotype is "enhanced" or
  "reduced").

  All operations on a \code{GeneSet} can be applied to a
  \code{GeneColorSet}; coloring can also be accesssed.
}
\section{Objects from the Class}{
  Construct a \code{GeneColorSet} with a \code{\link{GeneColorSet}}
  method. These methods are identical to those for \code{GeneSet},
  except they require an additional \code{phenotype} argument to
  specify the phenotype to which the genetic and phenotypic coloring
  apply. A \code{GeneColorSet} can be constructed from a \code{GeneSet}
  with \code{GeneColorSet(<GeneSet>, phenotype="<phenotype>").}
}
\section{Slots}{
  A \code{GeneColorSet} inherits all slots from
  \code{\linkS4class{GeneSet}}, and gains the following slots:
	 \describe{
    \item{\code{phenotype}:}{Object of class \code{"ScalarCharacter"}
      describing the phenotype for which this gene set is colored.}
    \item{\code{geneColor}:}{Object of class \code{"factor"} describing
      the coloring of each gene in the set. The lengths of
      \code{geneColor} and \code{gene} must be equal.}
    \item{\code{phenotypeColor}:}{Object of class \code{"factor"}
      describing the phenotypic coloring of each gene in the set. The
      lengths of \code{phenotypeColor} and \code{gene} must be equal.}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{GeneSet}"}, directly.
}
\section{Methods}{

  Methods unique to \code{GeneColorSet} include:
  \describe{
    \item{coloring}{\code{signature(object = "GeneColorSet")}: retrieve
      coloring as a \code{data.frame}. The row names of the data frame
      are the gene names; the columns are \code{geneColor} and
      \code{phenotypeColor}.}
    \item{coloring<-}{\code{signature(object = "GeneColorSet", value =
        "data.frame")}: use a \code{data frame} to assign coloring
      information. The \code{data.frame} must have the same number of
      rows as the \code{GeneColorSet} has genes (though see the exmaples
      below for flexible ways to alter coloring of a subset of
      genes). Row names of the \code{data.frame} correspond to gene
      names. The data frame has two columns, named \code{geneColor} and
      \code{phenotypeColor}. These must be of class \code{factor}.

      A typical use of \code{coloring<-} is to simultaneous extract,
      subset, and reassign the currenct coloring, e.g.,
      \code{coloring(<GeneColorSet>)[1:5,"geneColor"] <- "up"}; see the
      examples below.}
    \item{geneColor<-}{\code{signature(object = "GeneColorSet", value =
        "factor")}: assign gene colors. }
    \item{geneColor}{\code{signature(object = "GeneColorSet")}: retrieve
      gene colors as a \code{factor}.}
    \item{phenotypeColor<-}{\code{signature(object = "GeneColorSet",
        value = "factor")}: assign phenotype colors.}
    \item{phenotypeColor}{\code{signature(object = "GeneColorSet")}:
      retrieve phenotype colors as a \code{factor}.}
    \item{phenotype<-}{\code{signature(object = "GeneColorSet", value =
        "character")}: assign the phenotype from a single-element
      character vector.}
    \item{phenotype}{\code{signature(object = "GeneColorSet")}: retrieve
      the phenotype as a single-element \code{character}.}
  }

  \code{GeneColorSet} inherits all methods from class
  \code{GeneSet}. Methods with different behavior include
  \describe{
    \item{[}{
      \code{signature(x = "GeneSet", i="character")}
      \code{signature(x = "GeneSet", i="numeric")}: subset the gene set by
      index (\code{i="numeric"}) or gene value (\code{i="character"}). Genes
      are re-ordered as required. \code{geneColor} and \code{phenotypeColor} are
      subsetted as appropriate.}
    \item{[[}{\code{signature(x = "GeneSet")}: select a single gene from
      the gene set, returning a named character vector of \code{gene,
        geneColor, phenotypeColor}. Exact matches only.}
    \item{\$}{\code{signature(x = "GeneSet")}: select a single gene from
      the gene set, returning a named character vector of \code{gene,
        geneColor, phenotypeColor}. Provides partial matching into the
      list of genes.}
  }  

  Logical (set) operations \code{&, |, setdiff} warn if
  the \code{phenotype} \code{geneColor}, or \code{phenotypeColor}
  differs between sets; this implies coercion of factor levels, and the
  consequences should be carefully considered.
}
\author{Martin Morgan <mtmorgan@fhcrc.org>} 
\seealso{
  \code{\linkS4class{GeneSet}}.
}
\examples{
## Create a GeneColorSet from an ExpressionSet
data(sample.ExpressionSet)
gcs1 <- GeneColorSet(sample.ExpressionSet[100:109],
                     setName="color1", setIdentifier="200",
                     phenotype="imaginary")
gcs1
## or with color...
gcs2 <- GeneColorSet(sample.ExpressionSet[100:109],
                     setName="color1", setIdentifier="200",
                     phenotype="imaginary",
                     geneColor=factor(
                       rep(c("up", "down", "unchanged"),
                           length.out=10)),
                     phenotypeColor=factor(
                       rep(c("enhanced", "reduced"),
                           length.out=10)))
coloring(gcs2)

## recode geneColor of genes 1 and 4
coloring(gcs2)[c(1,4),"geneColor"] <- "down"
coloring(gcs2)
## reset, this time by gene name
coloring(gcs2)[c("31339_at", "31342_at"),"geneColor"] <- c("up", "up")
## usual 'factor' errors and warning apply:
coloring(gcs2)[c("31339_at", "31342_at"),"geneColor"] <- c("UP", "up")

gcs2[["31342_at"]]
try(gcs2[["31342_"]]) # no partial matching
gcs2$"31342" # 1 partial match ok
}
\keyword{classes}
