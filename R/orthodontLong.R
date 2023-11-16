#' Growth curve data on an orthodontic measurement
#'
#' The `orthodontLong` data.frame has 108 rows and 4 columns of the change in
#' an orthodontic measurement over time for several young subjects.  This
#' dataset was modified from the \code{\link[nlme]{Orthodont}} dataset in the
#' \code{\link[nlme]{nlme}} package.
#'
#' @format This `data.frame` contains the following columns:
#' \describe{
#'     \item{distance}{a numeric vector of distances from the pituitary to the
#'                     pterygomaxillary fissure (mm). These distances are
#'                     measured on x-ray images of the skull.}
#'     \item{age}{a factor indicating the ages of the subjects (yr).}
#'     \item{Subject}{a factor indicating the subject on which the measurement
#'                    was made. The levels are labelled `M01` to `M16` for the
#'                    males and `F01` to `F13` for the females.}
#'     \item{Sex}{a factor with levels `Male` and `Female`}
#'     }
#'
#' @details Investigators at the University of North Carolina Dental School
#'          followed the growth of 27 children (16 males, 11 females) from age 8
#'          until age 14. Every two years they measured the distance between the
#'          pituitary and the pterygomaxillary fissure, two points that are
#'          easily identified on x-ray exposures of the side of the head.
#'
#' @source {Potthoff, R. F. and Roy, S. N. (1964), "A generalized multivariate
#' analysis of variance model useful especially for growth curve problems",
#' Biometrika, 51, 313-326.}
#'
#' @examples
#' data(orthodontLong)
"orthodontLong"
