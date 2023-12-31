% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mmemmuris.R
\name{coefs}
\alias{coefs}
\title{Get Coefficients from Model}
\usage{
coefs(fMod)
}
\arguments{
\item{fMod}{A model of class \code{\link[nlme]{gls}}, \code{\link[nlme]{lme}},
or \code{mlm}.}
}
\value{
This function will return a matrix of parameter estimates for the model.
}
\description{
This function will return the (generalized) least squares estimates of
parameters in marginal and mixed-effects models and multivariate linear
models.
}
\details{
Marginal and mixed-effects model (fixed-effects) coefficients can be
calculated with generalized least squares as

\strong{\eqn{\hat{\beta}}} \eqn{ = (}\strong{X'V} \eqn{^ {-1}}\strong{X}\eqn{) ^ {-}}\strong{X'V} \eqn{^ {-1}}\strong{y}

where \strong{X} is the design matrix of the independent variables,
\strong{V} is the marginal variance-covariance matrix of the random/repeated
effects where \strong{V = ZGZ' + R} (\strong{Z} is the random effects design matrix,
\strong{G} is the corresponding (co)variances of the random effects matrix, and
\strong{R} represents the (co)variances of the repeated effects matrix), and \strong{y}
is a vector of observations.

Multivariate linear model coefficients can be calculated with ordinary least
squares as

\strong{\eqn{\hat{B}}} \eqn{ = (}\strong{X_b'X_b} \eqn{) ^ {-}}\strong{X_b'Y}

where \strong{X_b} is the design matrix of the between-subject effects and
\strong{Y} is a combined matrix of observations.

To get the multivariate linear model coefficients from a \code{\link[nlme]{gls}} or \code{\link[nlme]{lme}}
model, use the \link[mmemmuris:Hmatrix]{Hmatrix} function as it calls this
function for each level of the within-subject effect.
}
\examples{
# Marginal model
# Unstructured covariance
glsMod <- gls(Distance ~ Sex * Age,
              correlation =
                corSymm(form = ~ as.numeric(Age) | Subject),
              weights = varIdent(form = ~ 1 | Age),
              na.action = na.omit,
              data = orthodontLong)
# Model coefficients
coefs(glsMod)

# Linear mixed-effects model
# Unstructured covariance
lmeMod <- lme(Distance ~ Sex * Age,
              random = ~ Age | Subject,
              na.action = na.omit,
              data = orthodontLong)
# Fixed-effects model coefficients
coefs(lmeMod)

# Multivariate linear model
manovaMod <- manova(cbind(Distance8, Distance10, Distance12, Distance14) ~ Sex,
                    data = orthodontWide)
# Model coefficients
coefs(manovaMod)

# Get multivariate linear model coefficients from a univariate model
coefs(glsMod)[1:2, ]
orthodontLong$Age <- relevel(orthodontLong$Age, ref = "10")
coefs(update(glsMod, data = orthodontLong))[1:2, ]
orthodontLong$Age <- relevel(orthodontLong$Age, ref = "12")
coefs(update(glsMod, data = orthodontLong))[1:2, ]
orthodontLong$Age <- relevel(orthodontLong$Age, ref = "14")
coefs(update(glsMod, data = orthodontLong))[1:2, ]

}
\seealso{
{\code{\link[mmemmuris]{Hmatrix}}, \code{\link[mmemmuris]{Ematrix}}, \code{\link[mmemmuris]{waldF}}, \code{\link[mmemmuris]{uniTest}}}
}
