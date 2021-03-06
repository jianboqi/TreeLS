% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stem_points_methods.R
\name{stm.hough}
\alias{stm.hough}
\title{Stem denoising algorithm: Hough Transform}
\usage{
stm.hough(
  h_step = 0.5,
  max_d = 0.5,
  h_base = c(1, 2.5),
  pixel_size = 0.025,
  min_density = 0.1,
  min_votes = 3
)
}
\arguments{
\item{h_step}{\code{numeric} - height interval to perform point filtering/assignment/classification.}

\item{max_d}{\code{numeric} - largest tree diameter expected in the point cloud.}

\item{h_base}{\code{numeric} vector of length 2 - tree base height interval to initiate circle search.}

\item{pixel_size}{\code{numeric} - pixel side length to discretize the point cloud layers while performing the Hough Transform circle search.}

\item{min_density}{\code{numeric} - between 0 and 1 - minimum point density within a pixel evaluated on the Hough Transform - i.e. only \emph{dense} point clousters will undergo circle search.}

\item{min_votes}{\code{integer} - Hough Transform parameter - minimum number of circle intersections over a pixel to assign it as a circle center candidate.}
}
\description{
This function is meant to be used inside \code{\link{stemPoints}}. It applies an adapted version of the Hough Transform for circle search. Mode details are given in the sections below.
}
\section{\code{LAS@data} Special Fields}{


Meaninful new fields in the output:

\itemize{
\item \code{Stem}: \code{TRUE} for stem points
\item \code{Segment}: stem segment number (from bottom to top and nested with TreeID)
\item \code{Radius}: approximate radius of the point's stem segment estimated by the Hough Transform - always a multiple of the \code{pixel_size}
\item \code{Votes}: votes received by the stem segment's center through the Hough Transform
}
}

\section{Adapted Hough Transform}{


The Hough Transform circle search algorithm used in
TreeLS applies a constrained circle search on discretized 
point cloud layers. Tree-wise, the circle search is  
recursive, in which the search for circle parameters 
of a stem section is constrained to the 
\emph{feature space} of the stem section underneath it.
Initial estimates of the stem's \emph{feature space} 
are performed on a \emph{baselise} stem segment - i.e.
a low height interval where a tree's bole is expected  
to be clearly visible in the point cloud.
The algorithm is described in detail by Conto et al. (2017).

This adapted version of the algorithm is very robust against outliers, 
but not against forked or leaning stems.
}

\references{


Conto, T. et al., 2017. Performance of stem denoising and stem modelling algorithms on single tree point clouds from terrestrial laser scanning. Computers and Electronics in Agriculture, v. 143, p. 165-176.
}
