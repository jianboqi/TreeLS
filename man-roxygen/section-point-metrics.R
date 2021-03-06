#' @section List of available point metrics:
#' 
#' \loadmathjax
#' 
#' * \emph{EVi} = \emph{i}-th 3D eigen value 
#' 
#' * \emph{EV2Di} = \emph{i}-th 2D eigen value
#' 
#' \itemize{
#'    \item \code{N}: number of nearest neighbors
#'    \item \code{MinDist}: minimum distance among neighbors 
#'    \item \code{MaxDist}: maximum distance among neighbors
#'    \item \code{MeanDist}: mean distance
#'    \item \code{SdDist}: standard deviation of within neighborhood distances
#'    \item \code{Linearity}: linear saliency, \mjeqn{(EV_{1} + EV_{2}) / EV_{1}}{(EV1 + EV2) / EV1}
#'    \item \code{Planarity}: planar saliency, \mjeqn{(EV_{2} + EV_{3}) / EV_{1}}{(EV2 + EV3) / EV1}
#'    \item \code{Scattering}: \mjeqn{EV_{3} / EV_{1}}{EV3 / EV1}
#'    \item \code{Omnivariance}: \mjeqn{(EV_{2} + EV_{3}) / EV_{1}}{(EV2 + EV3) / EV1}
#'    \item \code{Anisotropy}: \mjeqn{(EV_{1} - EV_{3}) / EV_{1}}{(EV1 - EV3) / EV1}
#'    \item \code{Eigentropy}: \mjeqn{- \sum_{i=1}^{n=3} EV_{i} * ln(EV_{i})}{-sum(EV * ln(EV))}
#'    \item \code{EigenSum}: sum of eigenvalues, \mjeqn{\sum_{i=1}^{n=3} EV_{i}}{sum(EV)}
#'    \item \code{Curvature}: surface variation, \mjeqn{EV_{3} / EigenSum}{EV3 / EigenSum}
#'    \item \code{KnnRadius}: 3D neighborhood radius
#'    \item \code{KnnDensity}: 3D point density (N / sphere volume)
#'    \item \code{Verticality}: absolute vertical deviation, in degrees 
#'    \item \code{ZRange}: point neighborhood height difference
#'    \item \code{ZSd}: standard deviation of point neighborhood heights
#'    \item \code{KnnRadius2d}: 2D neighborhood radius 
#'    \item \code{KnnDensity2d}: 2D point density (N / circle area)
#'    \item \code{EigenSum2d}: sum of 2D eigenvalues, \mjeqn{\sum_{i=1}^{n=2} EV2D_{i}}{sum(EV2D)}
#'    \item \code{EigenRatio2d}: \mjeqn{EV2D_{2} / EV2D_{1}}{EV2D2 / EV2D1}
#'    \item \code{EigenValuei}: 3D eigenvalues
#'    \item \code{EigenVectorij}: 3D eigenvector coefficients, \emph{i}-th load of \emph{j}-th eigenvector
#' }
