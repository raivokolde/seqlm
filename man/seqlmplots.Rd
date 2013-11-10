\name{seqlmplots}
\alias{seqlmplots}
\title{Visualise the regions}
\usage{
  seqlmplots(segments, values, genome_information,
    annotation, expand = 100, ylim = extendrange(values),
    dir = NA, filetype = "png", ...)
}
\arguments{
  \item{segments}{selection of significant regions by
  \code{\link{seqlm}} function}

  \item{values}{same values matrix that was used in
  \code{seqlm}}

  \item{genome_information}{same genome_information object
  that was used in \code{seqlm}}

  \item{annotation}{same annotation vector that was used in
  \code{seqlm}}

  \item{expand}{number of basepairs to extend the region on
  plot}

  \item{ylim}{two element vector giving the lower and
  higher limit of the y axis}

  \item{dir}{directory where to put the images, of NA then
  plots are drawn into the plotting window}

  \item{filetype}{picture filetype}

  \item{...}{extra parameters to \code{\link{ggsave}}}
}
\description{
  Generate plots about the seqlm results
}
\details{
  The number of results from \code{\link{seqlm}} can be
  large and visualising all these regions might not be
  desirable. Therefore, it is advisable to filter the
  results befor plotting.
}
\author{
  Raivo Kolde <rkolde@gmail.com>
}

