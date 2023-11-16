#' Growth curve data on an orthodontic measurement
#'
#' The `orthodontWide` data.frame has 27 rows and 6 columns of the change in 
#' an orthodontic measurement over time for several young subjects.  This dataset
#' was modified from the \code{\link[nlme]{Orthodont}} dataset in the 
#' \code{\link[nlme]{nlme}} package.
#' 
#' @format This `data.frame` contains the following columns:
#' \describe{
#'     \item{Subject}{a factor indicating the subject on which the measurement 
#'                    was made. The levels are labelled `M01` to `M16` for the 
#'                    males and `F01` to `F13` for the females.}
#'     \item{Sex}{a factor with levels `Male` and `Female`}
#'     \item{distance8}{a numeric vector of distances from the pituitary to the 
#'                 pterygomaxillary fissure (mm) for age 8 subjects. These
#'                 distances are measured on x-ray images of the skull.}
#'     \item{distance10}{a numeric vector of distances from the pituitary to the 
#'                 pterygomaxillary fissure (mm) for age 10 subjects. These
#'                 distances are measured on x-ray images of the skull.}
#'     \item{distance12}{a numeric vector of distances from the pituitary to the 
#'                 pterygomaxillary fissure (mm) for age 12 subjects. These
#'                 distances are measured on x-ray images of the skull.}
#'     \item{distance14}{a numeric vector of distances from the pituitary to the 
#'                 pterygomaxillary fissure (mm) for age 14 subjects. These
#'                 distances are measured on x-ray images of the skull.}
#'     }
#'     
#' @details Investigators at the University of North Carolina Dental School 
#'          followed the growth of 27 children (16 males, 11 females) from age 8
#'          until age 14. Every two years they measured the distance between the
#'          pituitary and the pterygomaxillary fissure, two points that are 
#'          easily identified on x-ray exposures of the side of the head.
#'
#' @source {Potthoff, R. F. and Roy, S. N. (1964), "A generalized multivariate analysis of variance model useful especially for growth curve problems", Biometrika, 51, 313-326.}
#'
#' @examples
#' data(orthodontWide)
"orthodontWide"