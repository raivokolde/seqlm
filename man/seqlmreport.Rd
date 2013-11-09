\name{seqlmreport}
\alias{seqlmreport}
\title{Generate the HTML report for the seqlm results}
\usage{
  seqlmreport(segments, values, genome_information,
    annotation, ylim = extendrange(values), dir = NA,
    expand = 100, width = 8, height = 5, dpi = 100,
    main = "seqlm results")
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

  \item{dir}{directory where to put the page, if the
  directory does not exist it will be created}

  \item{width}{picture width in inches}

  \item{height}{picture height in inches}

  \item{dpi}{dots per inch, to calibrate the picture size
  in pixels}

  \item{main}{title for the raport}
}
\description{
  Generate the HTML report for the seqlm results
}
\author{
  Kaspar Martens <kmartens@ut.ee> Raivo Kolde
  <rkolde@gmail.com>
}

