#' Blood histamine measurements on dogs
#'
#' The `dogsLong` data.frame has 64 rows and 5 columns of log-histamine
#' measurements over time for several dogs.
#'
#' @format This `data.frame` contains the following columns:
#' \describe{
#'     \item{Dog}{a factor indicating the dog on which the measurement was made}
#'     \item{Drug}{a factor with levels `Morphine` and `Trimethaphan`}
#'     \item{Depleted}{a factor with levels `Y` and `N`}
#'     \item{Minutes}{a factor indicating the minutes after the drug is
#'     injected.}
#'     \item{logHistamine}{a numeric vector of log-transformed blood
#'     concentrations of histamine at 0, 1, 3, or 5 minutes after injection of
#'     the drug}
#'     }
#'
#' @details Sixteen dogs are randomly assigned to four groups.  Dogs in each
#'          group receive either morphine or trimethaphan (variable Drug) and
#'          have either depleted or intact histamine levels (variable Depleted)
#'          before receiving the drugs. The dependent variable is the log-blood
#'          concentration of histamine at 0, 1, 3, and 5 minutes after injection
#'          of the drug.
#'
#' @source {Cole, J. W. L., and Grizzle, J. E. (1966). "Applications of
#' Multivariate Analysis of Variance to Repeated Measures Experiments."
#' Biometrics 22:810-828.}
#' @source {\url{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.4/statug/statug_glm_examples07.htm}}
#'
#' @examples
#' # Marginal model results
#' # Create the grouping factor
#' dogsLong$Group <- factor(paste(dogsLong$Drug, dogsLong$Depleted))
#' glsModAA <- gls(logHistamine ~ Group * Minutes,
#'                 correlation =
#'                   corSymm(form = ~ as.numeric(Minutes) | Dog),
#'                 weights = varIdent(form = ~ 1 | Minutes),
#'                 na.action = na.omit,
#'                 data = dogsLong)
#' # McKeon F
#' hlTrace(glsModAA)
#' # Pillai-Samson F
#' hlTrace(glsModAA, approximation = "P")
#'
#' # Mixed-effects model results
#' # Old way
#' lmeModAA <- lme(logHistamine ~ Group * Minutes,
#'                 random = ~ Minutes | Dog,
#'                 weights = varIdent(form = ~ 1 | Minutes),
#'                 na.action = na.omit,
#'                 data = dogsLong)
#' # McKeon F
#' hlTrace(lmeModAA)
#' # Pillai-Samson F
#' hlTrace(lmeModAA, approximation = "P")
#'
"dogsLong"
