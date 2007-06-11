\name{xal}
\alias{xml}
\alias{getBroadSets}
\alias{toBroadXML}

\title{Read and Write XML Gene Sets}
\description{
  \code{getBroadSets} parses one or more XML files for gene sets. The
  file can reside locally or at a URL. The format followed is that defined
  by the Broad (below).

  \code{toBroadXML} creates Broad XML from \code{BroadCollection} gene
  sets.
}
\usage{
getBroadSets(file, ...)
toBroadXML(geneSet, con, ...)
}
\arguments{
  \item{file}{A file name or URL containing gene sets encoded following
    the Broad specification. Gene sets can be embedded in other
    elements, such as those for Broad's MSIGDB.}
  \item{geneSet}{A \code{\linkS4class{GeneSet}} with
    \code{collectionType} \code{BroadCollection} (to ensure that
    required information is available).}
  \item{con}{An optional file name or connection to receive output.}
  \item{\dots}{Further arguments passed to the underlying XML parser,
    particularly \code{file} used to specify an output \code{connection}
    for \code{toBroadXML}.}
}
\value{
  \code{getBroadSets} returns a \code{GeneSetCollection} of gene sets.

  \code{toBroadXML} returns a character vector of a single
  \code{GeneSet} or, if \code{file} is provided, writes the XML to a
  file.
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
## could be a url
fl <- system.file("extdata", "Broad.xml", package="GSEABase")
gss <- getBroadSets(fl) # GeneSetCollection of 2 sets
names(gss)
gss[[1]]

fl <- tempfile()
toBroadXML(gss[[1]], con=fl)
noquote(readLines(fl))
unlink(fl)

toBroadXML(gss[[1]]) # character vector
}
\keyword{manip}
