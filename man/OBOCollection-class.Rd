\name{OBOCollection-class}
\docType{class}
% OBOCollection-class
\alias{OBOCollection-class}
\alias{subsets}
\alias{subsets,OBOCollection-method}
\alias{[,OBOCollection,character,missing-method}
\alias{coerce,OBOCollection,graphNEL-method}
\alias{coerce,graphNEL,OBOCollection-method}
\alias{show,OBOCollection-method}

\title{Class "OBOCollection"}

\description{
  \code{OBOCollection} extends the \code{\link{GOCollection}} class, and
  is usually constructed from a file formated following the OBO file
  format. See \code{\linkS4class{CollectionType}} for general use of
  collections with gene sets.
}

\section{Objects from the Class}{
  Objects are instantiated with calls to \code{\link{OBOCollection}} or
  \code{\link{getOBOCollection}}.
}

\section{Slots}{
  \code{OBOCollection} extends \code{GOCollection} and
  \code{OBOCollection} has the following additional slots (these slots are
  NOT meant to be manipulated directly by the user):
  \describe{
    \item{\code{.stanza}:}{A \code{data.frame} representing the stanzas
      present in an OBO file. Row names of the data frame are
      unique stanza identifiers. The \code{value} column contains the
      stanza name (e.g., \sQuote{Term}, i.e., the stanza name associated
      with a GO identifier).}
    \item{\code{.subset}}{A \code{data.frame} representing (optional)
      subsets defined in the collection. Subsets are defined in the
      header of an OBO file with a subsetdef tag. Row names of the data
      frame are the subsetdef names; the \code{value} column contains
      the subset definition.}
    \item{\code{.kv}}{A \code{data.frame} representing key-value pairs
      in the OBO source file. The row names of the data frame correspond
      to lines in the OBO file. The \code{stanza_id} column indexes the
      row of \code{.stanza} describing the stanza in which the key-value
      pair occured. The remaining columns (\code{key}, \code{value})
      contain the parsed key and value.}
  }
}

\section{Methods}{
  \code{OBOCollection} has the following methods, in addition to those inherited from
  \code{\linkS4class{GOCollection}}.

  These methods list and select subsets of \code{OBOCollection}:
  \describe{
    \item{subsets}{\code{signature(object="OBOCollection",
        display="named")}: return a character vector of subsets present
      in \code{object}. Valid values for \code{display} are
      \sQuote{named} (a named character vector, with names equal to the
      names of the subsets and values the descriptions), \sQuote{full}
      (a character vector of name and description, with each pair
      formated into a single entry as \dQuote{name (description)}),
      \sQuote{key} (subset names), or \sQuote{value} (subset
      descriptions).}
    \item{[}{\code{signature(object="OBOCollection", i="character",
        j="missing", ...)}: return an
      \code{OBOCollection} by selecting just those subsets whose name
      matches the string(s) in \code{i}. This method calls the
      \code{[,GOCollection} method so, e.g., \code{evidenceCode} can be
      used to restricts which evidence codes the collection will
      identify.}
  }

  These methods coerce to and from \code{OBOCollection}:
  \describe{
    \item{as}{\code{signature(object="OBOCollection", "graphNEL")}:
      create a directed graph with nodes generated from
      \code{ids(object)} and edges from \code{is_a} relations of
      \code{object}.}
    \item{as}{\code{signature(object="graphNEL", "OBOCollection")}:
      create an \code{OBOCollection} with \code{ids} from the graph
      nodes, and edges from \code{inNodes(object)}.}
  }
}

\author{Martin Morgan <mtmorgan@fhcrc.org>}
\seealso{
  \code{\link{OBOCollection}} constructor; 
  \code{\linkS4class{CollectionType}} classes.
}
\references{
  \url{http://www.geneontology.org} for details of the OBO format.
}
\examples{
fl <- system.file("extdata", "goslim_plant.obo", package="GSEABase")
obo <- getOBOCollection(fl)
obo
subsets(obo)
obo["goslim_plant", evidenceCode="TAS"]
g <- as(obo["goslim_goa"], "graphNEL")
if (interactive() && require("Rgraphviz")) {
    plot(g)
}
}
\keyword{classes}
