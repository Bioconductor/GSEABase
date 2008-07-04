\name{import/export}
\alias{io}
\alias{getBroadSets}
\alias{asBroadUri}
\alias{toBroadXML}

\alias{toGmt}
\alias{getGmt}

\title{Read and write gene sets from Broad or GMT formats}
\description{

  \code{getBroadSets} parses one or more XML files for gene sets. The
  file can reside locally or at a URL. The format followed is that
  defined by the Broad (below).  \code{toBroadXML} creates Broad XML
  from \code{BroadCollection} gene sets.

  \code{toGmt} converts \code{GeneSetColletion} objects to a character
  vector representing the gene set collection in GMT
  format. \code{fromGmt} reads a file or other character vector into a
  \code{GeneSetColletion}.
}
\usage{
getBroadSets(uri, ...)
toBroadXML(geneSet, con, ...)
asBroadUri(name,
           base="http://www.broad.mit.edu/gsea/msigdb/cards")
getGmt(con, geneIdType=NullIdentifier(),
       collectionType=NullCollection(), sep="\t", ...)
toGmt(x, con, ...)
}
\arguments{
  \item{uri}{A file name or URL containing gene sets encoded following
    the Broad specification. For Broad sets, the uri can point to
    a MSIGDB.}
  \item{geneSet}{A \code{\linkS4class{GeneSet}} with
    \code{collectionType} \code{BroadCollection} (to ensure that
    required information is available).}
  \item{x}{A \code{\linkS4class{GeneSetCollection}} or other object for
    which a \code{toGmt} method is defined.}
  \item{con}{A (optional, in the case of \code{toXxx}) file name or
    connection to receive output.}
  \item{name}{A character vector of Broad gene set names, e.g.,
    \code{c('chr16q', 'GNF2_TNFSF10')}.}
  \item{base}{Base uri for finding Broad gene sets.}

  \item{geneIdType}{A constructor for the type of identifier the members
    of the gene sets represent. See \code{\link{GeneIdentifierType}} for
    more information.}

  \item{collectionType}{A constructor for the type of collection for the
    gene sets. See \code{\link{CollectionType}} for more information.}

  \item{sep}{The character string separating members of each gene set in
    the GMT file.}

  \item{\dots}{Further arguments passed to the underlying XML parser,
    particularly \code{file} used to specify an output \code{connection}
    for \code{toBroadXML}.}
}
\value{
  \code{getBroadSets} returns a \code{GeneSetCollection} of gene sets.

  \code{toBroadXML} returns a character vector of a single
  \code{GeneSet} or, if \code{con} is provided, writes the XML to a
  file.

  \code{asBroadUri} can be used to create URI names (to be used by
  \code{getBroadSets} of Broad files.

  \code{getGmt} returns a \code{GeneSetCollection} of gene sets.

  \code{toGmt} returns character vectors where each line represents a
  gene set. If \code{con} is provided, the result is written to the
  specified connection.

}

\references{
  \url{http://www.broad.mit.edu/gsea/}
}
\author{Martin Morgan <mtmrogan@fhcrc.org>}
\note{
  Actual Broad XML files differ from the DTD (e.g., an implied ','
  separator between genes in a set); we parse to and from files as they
  exists the actual files.
}
\seealso{
  \code{\linkS4class{GeneSetCollection}}
  \code{\linkS4class{GeneSet}}
}
\examples{
## 'fl' could also be a URI
fl <- system.file("extdata", "Broad.xml", package="GSEABase")
gss <- getBroadSets(fl) # GeneSetCollection of 2 sets
names(gss)
gss[[1]]

\dontrun{
## Download from the Broad
getBroadSets(asBroadUri(c('chr16q', 'GNF2_ZAP70')))
}

fl <- tempfile()
toBroadXML(gss[[1]], con=fl)
noquote(readLines(fl))
unlink(fl)

\dontrun{
toBroadXML(gss[[1]]) # character vector
}

fl <- tempfile()
toGmt(gss, fl)
getGmt(fl)
unlink(fl)
}
\keyword{manip}
