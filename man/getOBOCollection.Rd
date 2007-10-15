\name{getOBOCollection}
\alias{getOBOCollection}

\title{Read OBO-specified Gene Ontology Collections}
\description{
  \code{getOBOCollection} parses a uri (file or internet location)
  encoded following the OBO specification defined by the Gene Onotology
  consortium.
}
\usage{
getOBOCollection(uri, evidenceCode="ANY", ...)
}
\arguments{
  \item{uri}{A file name or URL containing gene sets encoded following
    the OBO specification.}
  \item{evidenceCode}{A character vector of evidence codes.}
  \item{\dots}{Further arguments passed to the
    \code{\link{OBOCollection}} constructor.}
}
\value{
  \code{getOBOCollection} returns an \code{OBOCollection} of gene
  sets. The gene set is constructed by parsing the file for \code{id}
  tags in \code{TERM} stanzas. The parser does not currently support all
  features of OBO, e.g., the ability to import additional files.
}

\references{
  \url{http://www.geneontology.org}
}
\author{Martin Morgan <mtmrogan@fhcrc.org>}
\seealso{
  \code{\linkS4class{OBOCollection}}, \code{\link{OBOCollection}}
}
\examples{
## 'fl' could also be a URI
fl <- system.file("extdata", "goslim_plant.obo", package="GSEABase")
getOBOCollection(fl) # GeneSetCollection of 2 sets

\dontrun{
## Download from the internet
fl <- "http://www.geneontology.org/GO_slims/goslim_plant.obo"
getOBOCollection(fl, evidenceCode="TAS")
}
}
\keyword{manip}
