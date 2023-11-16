#' Blood histamine measurements on dogs
#'
#' The `dogsWide` data.frame has 16 rows and 7 columns of log-histamine
#' measurements over time for several dogs.
#'
#' @format This `data.frame` contains the following columns:
#' \describe{
#'     \item{Dog}{a factor indicating the dog on which the measurement was made}
#'     \item{Drug}{a factor with levels `Morphine` and `Trimethaphan`}
#'     \item{Depleted}{a factor with levels `Y` and `N`}
#'     \item{logHistamine0}{a numeric vector of log-transformed blood
#'     concentrations of histamine at 0 minutes after injection of the drug}
#'     \item{logHistamine1}{a numeric vector of log-transformed blood
#'     concentrations of histamine at 1 minutes after injection of the drug}
#'     \item{logHistamine3}{a numeric vector of log-transformed blood
#'     concentrations of histamine at 3 minutes after injection of the drug}
#'     \item{logHistamine5}{a numeric vector of log-transformed blood
#'     concentrations of histamine at 5 minutes after injection of the drug}
#'     }
#'
#' @details Sixteen dogs are randomly assigned to four groups.  Dogs in each
#'          group receive either morphine or trimethaphan (variable Drug) and have
#'          either depleted or intact histamine levels (variable Depleted) before
#'          receiving the drugs. The dependent variable is the log-blood concentration
#'          of histamine at 0, 1, 3, and 5 minutes after injection of the drug.
#'
#' @source {Cole, J. W. L., and Grizzle, J. E. (1966). "Applications of
#' Multivariate Analysis of Variance to Repeated Measures Experiments."
#' Biometrics 22:810-828.}
#' @source {\url{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.4/statug/statug_glm_examples07.htm}}
#'
#' @examples
#' # MANOVA results
#' # Create the grouping factor
#' dogsWide$Group <- factor(paste(dogsWide$Drug, dogsWide$Depleted))
#'
#' # Fit the model
#' manovaMod <- manova(cbind(logHistamine0, logHistamine1, logHistamine3,
#'                           logHistamine5) ~ Group, data = na.omit(dogsWide))
#'
#' # Get multivariate tests
#' # McKeon F
#' hlTrace(manovaMod)
#' # Pillai-Samson F
#' hlTrace(manovaMod, approximation = "P")
"dogsWide"
