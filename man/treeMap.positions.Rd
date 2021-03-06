% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{treeMap.positions}
\alias{treeMap.positions}
\title{Convert a tree map to a 2D \code{data.table}}
\usage{
treeMap.positions(map, plot = TRUE)
}
\arguments{
\item{map}{object generated by \code{\link{treeMap}}.}

\item{plot}{\code{logical} - plot the tree map?}
}
\value{
signed \code{data.table} of tree IDs and XY coordinates.
}
\description{
Extracts the tree XY positions from a \emph{treeMap} output.
}
\examples{
file = system.file("extdata", "pine_plot.laz", package="TreeLS")
tls = readTLS(file) \%>\%
  tlsNormalize \%>\%
  tlsSample

x = plot(tls)

map = treeMap(tls, map.hough(h_step = 1, max_h = 4))
add_treeMap(x, map, color='red')

xymap = treeMap.positions(map)
}
