#' Terms Type
#'
#' @export
#' @description
#' This function will determine the between and within-subject effects from a
#' \code{\link[nlme]{gls}} or \code{\link[nlme]{lme}}
#' model object.
#' @param fMod A model of class \code{\link[nlme]{gls}} or \code{\link[nlme]{lme}}.
#' @returns \tabular{ll}{
#'    `between` \tab A character vector of between-subject effects. \cr
#'    \tab \cr
#'    `within` \tab A character vector of within-subject effects. \cr
#'    \tab \cr
#'    `Xbetween` \tab A sum-to-zero design matrix of between-subject effects.
#'    \cr \tab \cr
#' }
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#' termsType(glsMod)
#'
#' # Linear mixed-effects model
#' # Unstructured covariance
#' lmeMod <- lme(Distance ~ Sex * Age,
#'               random = ~ Age | Subject,
#'               na.action = na.omit,
#'               data = orthodontLong)
#' termsType(lmeMod)

termsType <- function(fMod){
  mmemmuris:::modelCheck(fMod)
  dataset <- mmemmuris:::getDataset(fMod)
  vars <- attr(terms(fMod), "term.labels")[attr(terms(fMod), "term.labels") %in% colnames(dataset)]
  groups <- paste(gsub("~", "", nlme::getGroupsFormula(fMod)), collapse = "")
  vars <- sapply(vars, function(x){
    nrow(na.omit(unique(dataset[, groups, drop = FALSE]))) == nrow(na.omit(unique(dataset[, c(groups, x)])))
  })
  # Classify variables as between or within-subject effects
  between <- names(vars)[vars]
  within <- names(vars)[vars == FALSE]

  if(length(within) == 0L)
    stop("There are no within-subject effects in the model.", call. = FALSE)
  if(length(within) > 1L)
    stop("Multiple within-subject effects are not supported.", call. = FALSE)
  dat <-
    unique(dataset[, c(groups, between), drop = FALSE])

  # Add the interactions
  modelTerms <- attr(terms(fMod), "term.labels")
  within <- modelTerms[sapply(strsplit(modelTerms, ":"), function(x) { any(x %in% within) })]
  between <- modelTerms[!(modelTerms %in% within)]

  formulae <- "~ 1"
  if(length(between) != 0L)
    formulae <- paste(formulae, "+", paste(between, collapse = "+"))

  Xbetween <- model.matrix(as.formula(formulae), data = dat)
  rownames(Xbetween) <- NULL
  if(length(between) > 0L){
    X <- model.matrix(as.formula(paste("~ ", paste(c(between, within), collapse = "+"))),
                      data = dataset)
  }else{
    X <- model.matrix(as.formula(paste("~ ", paste(within, collapse = "+"))),
                      data = dataset)
  }
  rownames(X) <- NULL
  return(list(between = between, within = within, X = X, Xbetween = Xbetween))
}

uniTest.mlm <- function(fMod, within = "Time",
                        epsilon = c("Greenhouse-Geisser", "Huynh-Feldt-Lecoutre")){
  epsilon <- match.arg(epsilon)
  Esingular <- FALSE
  Ematrices <- mmemmuris:::Ematrix.mlm(fMod)
  if(is.null(mmemmuris:::inverseMatrix(Ematrices$E)))
    Esingular <- TRUE
  if(Esingular)
    stop("The SSCP E matrix is singular.  Univariate tests not available.", call. = FALSE)

  unV <- Ematrices$EB / Ematrices$n

  dat <- na.omit(eval(fMod$call$data))
  formulae <- paste(attr(terms(fMod), "term.labels"), collapse = "+")
  if(formulae == "")
    formulae <- "1"
  dat$y <- rowSums(Ematrices$Y) / sqrt(ncol(Ematrices$Y))
  test <- aov(as.formula(paste("y ~", formulae)), data = dat)
  betweenTests <- as.data.frame(summary(test)[[1]])
  colnames(betweenTests) <- c("df", "Sum Sq", "Mean Sq", "F", "Pr(>F)")
  betweenTests <- cbind(betweenTests, mmemmuris::sigStars(betweenTests$`Pr(>F)`))
  colnames(betweenTests)[length(colnames(betweenTests))] <- betweenTests[length(rownames(betweenTests)), length(colnames(betweenTests))] <- ""

  rowNames <- gsub(" ", "", rownames(betweenTests), fixed = TRUE)

  dfN <- betweenTests$df[!(rowNames %in% "Residuals")]
  names(dfN) <- rowNames[!(rowNames %in% "Residuals")]

  ddf <- "between-within"
  type <- "1"
  betweenTests <- list(betweenTests = betweenTests, ddf = ddf, type = type)
  class(betweenTests) <- c("list", "betweenTests")

  Hmatrices <- mmemmuris:::Hmatrix.mlm(fMod, within = within)

  dfNum <- c(1, dfN) * nrow(Ematrices$E)
  dfDen <- rep(fMod$df.residual * nrow(Ematrices$E), length(Hmatrices$H))

  F <- ((Hmatrices$SSH / dfNum) / (Ematrices$SSE / dfDen))

  withinTests <- data.frame(SSH = Hmatrices$SSH, dfNum = dfNum, SSE = Ematrices$SSE, dfDen = dfDen,
                            F = F, p.value = 1 - pf(F, dfNum, dfDen))
  colnames(withinTests)[c(2, 4, 6)] <- c("Num df", "Den df", "Pr(>F)")
  withinTests <- cbind(withinTests, mmemmuris::sigStars(withinTests$`Pr(>F)`))
  colnames(withinTests)[length(colnames(withinTests))] <- ""

  sphericityTestsRes <- mmemmuris::sphericityTests(Ematrices)
  if((sphericityTestsRes$r - 1) == 1L){
    sphericityTests <- epsCorrect <- NULL
    correctMethod <- "Lower-Bound"
  }else{
    eps <- mmemmuris::epsilon(Ematrices, method = epsilon)
    if(epsilon == "Greenhouse-Geisser"){
      if(eps > 0.75)
        warning("The Huynh-Feldt-Lecoutre epsilon is preferred when the Greenhouse-Geisser epsilon is greater than 0.75.", call. = FALSE)
      epsCorrect <- data.frame(eps,
                               eps * as.numeric(withinTests$`Num df`),
                               eps * as.numeric(withinTests$`Den df`),
                               as.numeric(withinTests$F),
                               1 - pf(as.numeric(withinTests$F),
                                      eps * as.numeric(withinTests$`Num df`),
                                      eps * as.numeric(withinTests$`Den df`)))
      correctMethod <- "Greenhouse-Geisser"
    }else if(epsilon == "Huynh-Feldt-Lecoutre"){
      epsCorrect <- data.frame(eps,
                               min(eps, 1) * as.numeric(withinTests$`Num df`),
                               min(eps, 1) * as.numeric(withinTests$`Den df`),
                               as.numeric(withinTests$F),
                               1 - pf(as.numeric(withinTests$F),
                                      min(eps, 1) * as.numeric(withinTests$`Num df`),
                                      min(eps, 1) * as.numeric(withinTests$`Den df`)),
                               row.names = NULL)
      correctMethod <- "Huynh-Feldt-Lecoutre"
    }else{
      epsCorrect <- data.frame(1,
                               as.numeric(withinTests$`Num df`),
                               as.numeric(withinTests$`Den df`),
                               as.numeric(withinTests$F),
                               1 - pf(as.numeric(withinTests$F),
                                      as.numeric(withinTests$`Num df`),
                                      as.numeric(withinTests$`Den df`)),
                               row.names = NULL)
      correctMethod <- "Lower-Bound"
    }
  }

  withinTests <- list(withinTests = withinTests, ddf = ddf, type = type)
  class(withinTests) <- c("list", "withinTests")

  if((sphericityTestsRes$r - 1) != 1L){
    rownames(epsCorrect) <- rownames(withinTests$withinTests)
    colnames(epsCorrect) <- c("epsilon", "Num df", "Den df", "F", "Pr(>F)")
    epsCorrect <- cbind(epsCorrect, mmemmuris::sigStars(epsCorrect$`Pr(>F)`))
    colnames(epsCorrect)[length(colnames(epsCorrect))] <- ""
    epsCorrect <- list(epsCorrect = epsCorrect, correctMethod = correctMethod)
    class(epsCorrect) <- c("list", "epsCorrect")
  }

  csV <- mmemmuris:::Vmatrix.mlm(fMod, within = within)$csV
  csV <- list(V = csV)
  class(csV) <- c("list", "Vmatrix.mlm")

  uniTest <- list(H = Hmatrices$H, E = Ematrices$E,
                  betweenTests = betweenTests,
                  sphericityTests = sphericityTestsRes,
                  withinTests = withinTests,
                  epsCorrect = epsCorrect,
                  unV = unV, csV = csV, n = nrow(dat))
  class(uniTest) <- c("list", "uniTest", "mlm")
  return(uniTest)
}

uniTest.ulm <- function(fMod, individual = NULL, epsilon = c("Greenhouse-Geisser", "Huynh-Feldt-Lecoutre"), refit = TRUE, ss = NULL){
  epsilon <- match.arg(epsilon)
  ddf <- "between-within"
  type <- "1"
  Esingular <- FALSE
  if(mmemmuris::covStruct(fMod) %in% c("cs", "other"))
    stop("Only unstructured covariance models are allowed.", call. = FALSE)
  tt <- mmemmuris::termsType(fMod)
  data <- mmemmuris:::getDataset(fMod)
  Ematrices <- mmemmuris:::Ematrix.ulm(fMod, individual, ss)
  unV <- Ematrices$Vmatrix
  E <- Ematrices$E
  if(is.null(mmemmuris:::inverseMatrix(E)))
    Esingular <- TRUE
  M <- Ematrices$M
  n <- Ematrices$n

  if(refit == TRUE | Esingular == TRUE){
    # Update the model to compound symmetry covariance
    if(class(fMod) == "gls"){
      fMod <- update(fMod,
                     correlation = nlme::corCompSymm(form = as.formula(paste("~ 1 | ", as.character(nlme::getGroupsFormula(fMod))[2]))),
                     weights = NULL,
                     na.action = na.omit,
                     method = "REML")

    }else{
      fMod <- update(fMod,
                     random = as.formula(paste("~ 1 | ", as.character(nlme::getGroupsFormula(fMod))[2], "")),
                     correlation = NULL,
                     weights = NULL,
                     na.action = na.omit,
                     method = "REML")
    }
    L <- mmemmuris:::Lmatrix.ulm(fMod)
    wald <- mmemmuris::waldF(fMod, L)
    wald$DenDF <- mmemmuris::ddfBW(fMod)$`Den df`
    wald$p.value <- 1 - pf(wald$F, wald$NumDF, wald$DenDF)
    wald <- as.data.frame(do.call("cbind", wald))
    waldB <- wald[tt$between,, drop = FALSE]
    waldW <- wald[tt$within,, drop = FALSE]
    waldW$SSE <- Ematrices$SSE
    waldW$SSH <- (waldW$F * waldW$NumDF * waldW$SSE) / waldW$DenDF
    colnames(waldB)[c(2, 4:5)] <- colnames(waldW)[c(2, 4:5)] <- c("Num df", "Den df", "Pr(>F)")
    typeIB <- waldB[, c(2, 4, 1, 5)]
    typeIW <- waldW[, c(7, 2, 6, 4, 1, 5)]
    # Calculate Greenhouse-Geisser and Huynh-Feldt(-Lecoutre) epsilons for
    # violations of sphericity
    typeIB <- cbind(typeIB, mmemmuris::sigStars(typeIB$`Pr(>F)`))
    colnames(typeIB)[length(colnames(typeIB))] <- typeIB[length(rownames(typeIB)), length(colnames(typeIB))] <- ""

    typeIB <- list(betweenTests = typeIB, ddf = ddf, type = type)

    if(Esingular == TRUE){
      sphericityTests <- epsCorrect <- NULL
      correctMethod <- ""
    }else{
      sphericityTests <- mmemmuris::sphericityTests(Ematrices)
      if((sphericityTests$r - 1) == 1L){
        sphericityTests <- epsCorrect <- NULL
        correctMethod <- "Lower-Bound"
      }else{
        eps <- mmemmuris::epsilon(Ematrices, method = epsilon)
        if(epsilon == "Greenhouse-Geisser"){
          epsCorrect <- data.frame(eps,
                                   eps * as.numeric(typeIW$`Num df`),
                                   eps * as.numeric(typeIW$`Den df`),
                                   as.numeric(typeIW$F),
                                   1 - pf(as.numeric(typeIW$F),
                                          eps * as.numeric(typeIW$`Num df`),
                                          eps * as.numeric(typeIW$`Den df`)))
          rownames(epsCorrect) <- rownames(typeIW)
          colnames(epsCorrect) <- c("epsilon", "Num df", "Den df", "F", "Pr(>F)")
          epsCorrect <- cbind(epsCorrect, mmemmuris::sigStars(epsCorrect$`Pr(>F)`))
          colnames(epsCorrect)[length(colnames(epsCorrect))] <- ""
          correctMethod <- "Greenhouse-Geisser"
        }else if(epsilon == "Huynh-Feldt-Lecoutre"){
          epsCorrect <- data.frame(eps,
                                   min(eps, 1) * as.numeric(typeIW$`Num df`),
                                   min(eps, 1) * as.numeric(typeIW$`Den df`),
                                   as.numeric(typeIW$F),
                                   1 - pf(as.numeric(typeIW$F),
                                          min(eps, 1) * as.numeric(typeIW$`Num df`),
                                          min(eps, 1) * as.numeric(typeIW$`Den df`)),
                                   row.names = NULL)
          rownames(epsCorrect) <- rownames(typeIW)
          colnames(epsCorrect) <- c("epsilon", "Num df", "Den df", "F", "Pr(>F)")
          epsCorrect <- cbind(epsCorrect, mmemmuris::sigStars(epsCorrect$`Pr(>F)`))
          colnames(epsCorrect)[length(colnames(epsCorrect))] <- ""
          correctMethod <- "Huynh-Feldt-Lecoutre"
        }else{
          epsCorrect <- data.frame(1,
                                   as.numeric(typeIW$`Num df`),
                                   as.numeric(typeIW$`Den df`),
                                   as.numeric(typeIW$F),
                                   1 - pf(as.numeric(typeIW$F),
                                          as.numeric(typeIW$`Num df`),
                                          as.numeric(typeIW$`Den df`)),
                                   row.names = NULL)
          rownames(epsCorrect) <- rownames(typeIW)
          colnames(epsCorrect) <- c("epsilon", "Num df", "Den df", "F", "Pr(>F)")
          epsCorrect <- cbind(epsCorrect, mmemmuris::sigStars(epsCorrect$`Pr(>F)`))
          colnames(epsCorrect)[length(colnames(epsCorrect))] <- ""
          correctMethod <- "Lower-Bound"
        }
        epsCorrect <- list(epsCorrect = epsCorrect, correctMethod = correctMethod)
        class(epsCorrect) <- c("list", "epsCorrect")
      }
    }

      typeIW <- cbind(typeIW, mmemmuris::sigStars(typeIW$`Pr(>F)`))
      colnames(typeIW)[length(colnames(typeIW))] <- ""

      typeIW <- list(withinTests = typeIW, ddf = ddf, type = type)
      class(typeIW) <- c("list", "withinTests")

    Lwithin <- lapply(mmemmuris:::Lmatrix.ulm(fMod)[tt$within], t)
    csV <- mmemmuris:::Vmatrix.ulm(fMod, individual = individual)
    coefs <- mmemmuris::coefs(fMod)
    X <- model.matrix(fMod, data = mmemmuris:::getDataset(fMod))
    uniTest <- list(E = E, SSE = Ematrices$SSE,
                    betweenTests = typeIB,
                    sphericityTests = sphericityTests,
                    withinTests = typeIW,
                    epsCorrect = epsCorrect,
                    unV = unV, M = M, n = n,
                    Lwithin = Lwithin, csV = csV, coefs = coefs, X = X)
    class(uniTest) <- c("list", "uniTest", "ulm")
    return(uniTest)
  }else{
    Y <- do.call("rbind", tapply(mmemmuris:::getDataset(fMod)[[attr(getResponse(fMod), "label")]], mmemmuris:::getDataset(fMod)[[gsub("~", "", deparse(getGroupsFormula(fMod)))]], function(x) { x }))
    Y <- Y[match(unique(nlme::getGroups(fMod)), rownames(Y)), ]
    dat <- as.data.frame(unique(mmemmuris:::getDataset(fMod)[, c(gsub("~", "", deparse(getGroupsFormula(fMod))), tt$between[tt$between %in% colnames(mmemmuris:::getDataset(fMod))])]))
    dat$y <- rowSums(Y) / sqrt(ncol(Y))
    formulae <- paste(attr(terms(fMod), "term.labels")[attr(terms(fMod), "term.labels") %in% tt$between], collapse = "+")
    if(formulae == "")
      formulae <- "1"
    test <- aov(as.formula(paste("y ~", formulae)), data = dat)
    typeIB <- as.data.frame(summary(test)[[1]])
    colnames(typeIB)[c(1, 4)] <- c("df", "F")
    typeIB <- cbind(typeIB, mmemmuris::sigStars(typeIB$`Pr(>F)`))
    colnames(typeIB)[length(colnames(typeIB))] <- typeIB[length(rownames(typeIB)), length(colnames(typeIB))] <- ""
    typeIB <- list(betweenTests = typeIB, ddf = "between-within", type = type)
    class(typeIB) <- c("list", "betweenTests")

    H <- mmemmuris::Hmatrix(fMod)
    ddf <- mmemmuris::ddfBW(fMod)
    dfNum <- sapply(mmemmuris:::Lmatrix.ulm(fMod), Matrix::rankMatrix)[ddf$within]
    dfDen <- ddf$rDenDF - ddf$bDenDF
    Fstat <- ((H$SSH / dfNum) / (Ematrices$SSE / dfDen))
    typeIW <- data.frame(SSH = H$SSH, Numdf = dfNum, SSE = Ematrices$SSE, Dendf = dfDen, F = Fstat, p.value = 1 - pf(Fstat, dfNum, dfDen))
    colnames(typeIW)[c(2, 4, 6)] <- c("Num df", "Den df", "Pr(>F)")

    typeIW <- cbind(typeIW, mmemmuris::sigStars(typeIW$`Pr(>F)`))
    colnames(typeIW)[length(colnames(typeIW))] <- ""
    rownames(typeIW) <- tt$within

    typeIW <- list(withinTests = typeIW, ddf = "between-within", type = type)
    class(typeIW) <- c("list", "withinTests")

    if(Esingular == TRUE){
      sphericityTests <- epsCorrect <- NULL
      correctMethod <- ""
    }else{
      sphericityTests <- mmemmuris::sphericityTests(Ematrices)
      if((sphericityTests$r - 1) == 1L){
        sphericityTests <- epsCorrect <- NULL
        correctMethod <- "Lower-Bound"
      }else{
        eps <- mmemmuris::epsilon(Ematrices, method = epsilon)
        if(epsilon == "Greenhouse-Geisser"){
          epsCorrect <- data.frame(eps,
                                   eps * as.numeric(typeIW$withinTests$`Num df`),
                                   eps * as.numeric(typeIW$withinTests$`Den df`),
                                   as.numeric(typeIW$withinTests$F),
                                   1 - pf(as.numeric(typeIW$withinTests$F),
                                          eps * as.numeric(typeIW$withinTests$`Num df`),
                                          eps * as.numeric(typeIW$withinTests$`Den df`)))
          rownames(epsCorrect) <- rownames(typeIW$withinTests)
          colnames(epsCorrect) <- c("epsilon", "Num df", "Den df", "F", "Pr(>F)")
          epsCorrect <- cbind(epsCorrect, mmemmuris::sigStars(epsCorrect$`Pr(>F)`))
          colnames(epsCorrect)[length(colnames(epsCorrect))] <- ""
          rownames(epsCorrect) <- tt$within
          correctMethod <- "Greenhouse-Geisser"
        }else if(epsilon == "Huynh-Feldt-Lecoutre"){
          epsCorrect <- data.frame(eps,
                                   min(eps, 1) * as.numeric(typeIW$withinTests$`Num df`),
                                   min(eps, 1) * as.numeric(typeIW$withinTests$`Den df`),
                                   as.numeric(typeIW$withinTests$F),
                                   1 - pf(as.numeric(typeIW$withinTests$F),
                                          min(eps, 1) * as.numeric(typeIW$withinTests$`Num df`),
                                          min(eps, 1) * as.numeric(typeIW$withinTests$`Den df`)),
                                   row.names = NULL)
          rownames(epsCorrect) <- rownames(typeIW$withinTests)
          colnames(epsCorrect) <- c("epsilon", "Num df", "Den df", "F", "Pr(>F)")
          epsCorrect <- cbind(epsCorrect, mmemmuris::sigStars(epsCorrect$`Pr(>F)`))
          colnames(epsCorrect)[length(colnames(epsCorrect))] <- ""
          rownames(epsCorrect) <- tt$within
          correctMethod <- "Huynh-Feldt-Lecoutre"
        }else{
          epsCorrect <- data.frame(1,
                                   as.numeric(typeIW$withinTests$`Num df`),
                                   as.numeric(typeIW$withinTests$`Den df`),
                                   as.numeric(typeIW$withinTests$F),
                                   1 - pf(as.numeric(typeIW$withinTests$F),
                                          as.numeric(typeIW$withinTests$`Num df`),
                                          as.numeric(typeIW$withinTests$`Den df`)),
                                   row.names = NULL)
          rownames(epsCorrect) <- rownames(typeIW$withinTests)
          colnames(epsCorrect) <- c("epsilon", "Num df", "Den df", "F", "Pr(>F)")
          epsCorrect <- cbind(epsCorrect, mmemmuris::sigStars(epsCorrect$`Pr(>F)`))
          colnames(epsCorrect)[length(colnames(epsCorrect))] <- ""
          rownames(epsCorrect) <- tt$within
          correctMethod <- "Lower-Bound"
        }
        epsCorrect <- list(epsCorrect = epsCorrect, correctMethod = correctMethod)
        class(epsCorrect) <- c("list", "epsCorrect")
      }
    }

    Lwithin <- lapply(mmemmuris:::Lmatrix.ulm(fMod)[tt$within], t)
    coefs <- mmemmuris::coefs(fMod)
    X <- model.matrix(fMod, data = mmemmuris:::getDataset(fMod))
    uniTest <- list(E = E, SSE = Ematrices$SSE,
                    betweenTests = typeIB,
                    sphericityTests = sphericityTests,
                    withinTests = typeIW,
                    epsCorrect = epsCorrect,
                    unV = unV, M = M, n = n,
                    Lwithin = Lwithin, coefs = coefs, X = X)
    class(uniTest) <- c("list", "uniTest", "ulm")
    return(uniTest)
  }
}

#' Univariate Tests for RM ANOVA, Mixed-Effects Models, or Marginal Models
#'
#' @export
#'
#' @description
#' This function will calculate univariate tests for repeated measures ANOVA
#' (\code{mlm}), marginal (\code{\link[nlme]{gls}}) and mixed-effects
#' (\code{\link[nlme]{lme}}) models.  This includes Mauchly's test of sphericity
#' for within-subject effects and epsilon corrections for a departure from
#' sphericity.
#'
#' @details
#'
#' **H and E matrices**
#'
#' In RM ANOVA, there are within-subject effects and (optionally) between-subject
#' effects.  An ANOVA can be done on a transformed dependent variable denoted as
#'
#' \deqn{\Sigma(Y_i) / \sqrt r  for i = 1, \dots, r}
#'
#' where the numerator is the sum of the repeated measures responses for each
#' subject and the denominator is the square root of the number of repeated
#' measurements (\eqn{r}). The between-subject effects are the independent variables in the ANOVA to
#' obtain the between-subject effects tests.
#'
#' For within-subject effects, the SSCP error matrix is
#'
#' **E = M'** \eqn{(}**Y'Y**\eqn{-} **\eqn{\hat{B}}'** \eqn{(}**X'X**\eqn{)}**\eqn{\hat{B}}** \eqn{)}**M**
#'
#' and the SSCP hypothesis matrix is
#'
#' **H** \eqn{=} **M'** \eqn{(}**L\eqn{\hat{B}}** \eqn{)}'\eqn{(}**L**\eqn{(}**X'X**\eqn{) ^ {-}}**L'** \eqn{) ^ {-1}}**L\eqn{\hat{B}}M**
#'
#' where **M** is a sum-to-zero contrast matrix, **Y** is a combined
#' matrix of observations, **\eqn{\hat{B}}** is a matrix of combined coefficients from the
#' ordinary least squares models, and **X** is the design matrix of the
#' independent variables.
#'
#' For marginal and mixed effects models, the **E** matrix is calculated from
#' the marginal variance-covariance matrix of the random/repeated effects
#' where **V = ZGZ' + R** (**Z** is the random effects design matrix, **G** is
#' the corresponding (co)variances of the random effects matrix, and **R**
#' represents the (co)variances of the repeated effects matrix) as
#'
#' **E =** \eqn{n}**M'VM**
#'
#' and the **H** matrix is calculated with the same formula as above.  The **\eqn{\hat{B}}**
#' matrix is calculated by switching the reference levels of the within-subject
#' effect, refitting the models, and extracting the model coefficients.  The
#' contrast coefficients matrix (**L**) is calculated from the **LU**
#' decomposition of the crossproducts matrix (**X'X**).
#'
#' The F-statistics for the within-subject terms can be
#' calculated as \deqn{F = (SS_H / df_{Num}) / (SS_E / df_{Den})}
#'
#' where \eqn{SS_E = trace(}**E**\eqn{(}**M'M**\eqn{) ^ {-1})} is the error sum of
#' squares, the sum of squares for the hypothesis being tested is
#' \eqn{SS_H = trace(}**H**\eqn{(}**M'M**\eqn{) ^ {-1})}, and the numerator and
#' denominator degrees of freedom are denoted as \eqn{df_{Num}} and
#' \eqn{df_{Den}}, respectively.
#'
#' **Univariate Approach**
#'
#' The marginal or mixed-effects model is refit with a compound symmetry
#' covariance structure for the univariate tests with the `refit = TRUE`
#' argument.
#'
#' A Wald-type quadratic form is used to test each of the model terms.  The
#' Wald-type quadratic form for the F-tests is
#'
#' \eqn{F = (}**\eqn{\hat{\beta}}'L**\eqn{(}**L'** \eqn{(}**X'V**\eqn{ ^ {-1}}**X**\eqn{) ^ {-}}**L**\eqn{) ^ {-1}}**L'\eqn{\hat{\beta}}** \eqn{) / rank(}**L**\eqn{)}
#'
#' where **\eqn{\hat{\beta}}** is a vector of fixed effects coefficients, **L** is
#' a contrast coefficients matrix, **X** is the fixed effects design matrix, and
#' **V** is the marginal variance-covariance matrix
#' where **V = ZGZ' + R** (fit with restricted maximum likelihood).  **Z** is
#' the random effects design matrix, **G** is the corresponding (co)variances of
#' the random effects matrix, and **R** represents the (co)variances of the
#' repeated effects matrix.
#'
#' The F-statistics for the terms in the **H** and **E** matrices formula above
#' can be calculated as \deqn{F = (SS_H / df_{Num}) / (SS_E / df_{Den})}
#' With this information, we can backsolve for the pseudo-sum of squares for the
#' hypothesis matrices in marginal and mixed-effects models as
#' \deqn{SS_H = (F(df_{Num})(SS_E)) / df_{Den}} to reproduce this result in the
#' univariate tests table.
#'
#' **Epsilon Corrections and Mauchly's Sphericity Test**
#'
#' The within-subject SSCP **E** matrix is used to calculate the
#' Greenhouse-Geisser and Huynh-Feldt-Lecoutre epsilons for univariate
#' adjustments for violations of sphericity as well as the sphericity test
#' itself.
#'
#' The Greenhouse-Geisser epsilon is calculated as
#'
#' \eqn{\epsilon_{GG} = ((\Sigma\lambda_i of} **E**\eqn{(}**M'M**\eqn{) ^ {-1}) ^ 2) / (r - 1)(\Sigma\lambda_i ^ 2 of} **E**\eqn{(}**M'M**\eqn{) ^ {-1})}
#'
#' The Huynh-Feldt-Lecoutre epsilon is calculated as
#'
#' \eqn{\epsilon_{HFL} = ((v + 1)(r - 1)\epsilon_{GG} - 2) / ((r - 1)(v - (r - 1)\epsilon_{GG}))}
#'
#' Mauchly's test for sphericity is calculated for within-subject effects as
#'
#' \eqn{W = (\Pi\lambda_i of} **E**\eqn{(}**M'M**\eqn{) ^ {-1}) / (((\Sigma\lambda_i of} **E**\eqn{(}**M'M**\eqn{) ^ {-1}) / (r - 1)) ^ (r - 1))}
#'
#' Let \eqn{r} be equal to the number of repeated measurements, \eqn{v}
#' be the residual degrees of freedom (between-subject design matrix), and the correction factor equal
#' \deqn{\rho = 1 - (2(r - 1) ^ 2 + (r - 1) + 2) / (6(r - 1)v)}  The \eqn{\chi ^ 2} statistic is
#' equal to the product \eqn{-\rho v log(W)} with degrees of freedom
#' \eqn{(r(r - 1) / 2) - 1}.
#'
#' @param fMod A model of class \code{\link[nlme]{gls}}, \code{\link[nlme]{lme}}, or \code{mlm}.
#' @param individual The subject whose marginal variance-covariance matrix
#' should be used, where **V = ZGZ' + R** (**Z** is the random effects design
#' matrix, **G** is the corresponding (co)variances of the random effects matrix,
#' and **R** represents the (co)variances of the repeated effects matrix).  The
#' default is the first subject in the factor.
#' @param within The name that should be used for the within-subject effect.
#' The default is "Time".
#' @param epsilon Epsilon degrees of freedom adjustment used.  Options are
#' "Greenhouse-Geisser" and "Huynh-Feldt-Lecoutre".
#' @param refit The marginal or mixed-effects model is refit with a compound symmetry
#' covariance structure instead of calculating the statistics with the **H** and
#' **E** matrices.
#'
#' @references {\url{https://www.lesahoffman.com/PSYC943/mv12psyc943_lecture13.pdf}}
#' @references {\url{https://go.documentation.sas.com/doc/en/pgmsascdc/9.4_3.3/statug/statug_glm_details46.htm}}
#' @references {\url{https://support.sas.com/rnd/app/stat/papers/mixedglm.pdf}}
#' @references {\url{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.4/statug/statug_introreg_sect038.htm}}
#'
#' @returns \tabular{ll}{
#'    `EB` \tab A sum of squares cross products (SSCP) error matrix for
#'    between-subject effects. \cr
#'    \tab \cr
#'    `E` \tab A sum of squares cross products (SSCP) error matrix for
#'    within-subject effects. \cr
#'    \tab \cr
#'    `SSE` \tab The error sum of squares for within-subject effects calculated
#'    as \eqn{SS_E = trace(}**E**\eqn{(}**M'M**\eqn{) ^ {-1})}. \cr
#'    \tab \cr
#'    `spher` \tab A `data.frame` containing Mauchly's test for sphericity
#'    results. \cr
#'    \tab \cr
#'    `ggEpsilon` \tab The Greenhouse-Geisser epsilon value. \cr
#'    \tab \cr
#'    `hflEpsilon` \tab The Huynh-Feldt-Lecoutre epsilon value. \cr
#'    \tab \cr
#'    `V` \tab The **V** matrix from the marginal or mixed-effects model. \cr
#'    \tab \cr
#'    `M` \tab A sum-to-zero contrast matrix. \cr
#'    \tab \cr
#'    `n` \tab The sample size. \cr
#'    \tab \cr
#' }
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#' uniTest(glsMod)
#'
#' # Linear mixed-effects model
#' # Unstructured covariance
#' lmeMod <- lme(Distance ~ Sex * Age,
#'               random = ~ Age | Subject,
#'               na.action = na.omit,
#'               data = orthodontLong)
#' uniTest(lmeMod)
#'
#' # Multivariate linear model
#' manovaMod <-
#'   manova(cbind(Distance8, Distance10, Distance12, Distance14) ~ Sex,
#'          data = orthodontWide)
#' uniTest(manovaMod)
#'
#' @seealso
#' {\code{\link[mmemmuris]{Hmatrix}}, \code{\link[mmemmuris]{coefs}}, \code{\link[mmemmuris]{Vmatrix}}, \code{\link[mmemmuris]{Ematrix}}, \code{\link[mmemmuris]{waldF}}, \code{\link[mmemmuris]{sphericityTests}}, \code{\link[mmemmuris]{epsilon}}, \code{\link[mmemmuris]{ddfResidual}}, \code{\link[mmemmuris]{ddfBW}}}

uniTest <- function(fMod, individual = NULL, within = "Time", epsilon = c("Greenhouse-Geisser", "Huynh-Feldt-Lecoutre"), refit = TRUE, ss = NULL){
  if(any(class(fMod) %in% c("gls", "lme"))){
    epsilon <- match.arg(epsilon)
    uniTests <- mmemmuris:::uniTest.ulm(fMod, individual, epsilon, refit, ss)
  }else if(any(class(fMod) %in% "mlm")){
    epsilon <- match.arg(epsilon)
    uniTests <- mmemmuris:::uniTest.mlm(fMod, within, epsilon)
  }else
    stop('Please provide a model of class "gls", "lme", or "mlm".', call. = FALSE)
  return(uniTests)
}

wilksLambda.ulm <- function(fMod, individual = NULL, approximation = c("Rao", "Bartlett", "LR"), LR = TRUE, ss = NULL){
  approximation <- match.arg(approximation)
  tt <- mmemmuris::termsType(fMod)
  wilks <- mmemmuris:::lambda.ulm(fMod, individual, LR = LR, ss = ss)

  if(approximation == "Rao"){
    raoFApprox <- mmemmuris:::raoF(wilks)

    withinTests <- raoFApprox$wilksWithin

    if(length(raoFApprox$wilksBetween$wilks) > 0L){
      betweenTests <- raoFApprox$wilksBetween
      wilks <- list(betweenTests = betweenTests,
                    withinTests = withinTests,
                    LR = wilks$LR)
      class(wilks) <- c("list", "multiTest", "wilks", "ulm")
    }else{
      wilks <- list(withinTests = withinTests,
                    LR = wilks$LR)
      class(wilks) <- c("list", "multiTest", "wilks", "ulm")
    }
    return(wilks)
  }else if(approximation == "LR"){
    lrChiApprox <- mmemmuris:::lrChi(wilks)

    withinTests <- lrChiApprox$wilksWithin

    if(length(lrChiApprox$wilksBetween$wilks) > 0L){
      betweenTests <- lrChiApprox$wilksBetween
      wilks <- list(betweenTests = betweenTests,
                    withinTests = withinTests,
                    LR = wilks$LR)
      class(wilks) <- c("list", "multiTest", "wilks", "ulm")
    }else{
      wilks <- list(withinTests = withinTests,
                    LR = wilks$LR)
      class(wilks) <- c("list", "multiTest", "wilks", "ulm")
    }
    return(wilks)
  }else if(approximation == "Bartlett"){
    BChiApprox <- mmemmuris:::BChi(wilks)

    withinTests <- BChiApprox$wilksWithin

    if(length(BChiApprox$wilksBetween$wilks) > 0L){
      betweenTests <- BChiApprox$wilksBetween
      wilks <- list(betweenTests = betweenTests,
                    withinTests = withinTests,
                    LR = wilks$LR)
      class(wilks) <- c("list", "multiTest", "wilks", "ulm")
    }else{
      wilks <- list(withinTests = withinTests,
                    LR = wilks$LR)
      class(wilks) <- c("list", "multiTest", "wilks", "ulm")
    }
    return(wilks)
  }
}

wilksLambda.mlm <- function(fMod, approximation = c("Rao", "Bartlett", "LR"), within = "Time"){
  approximation <- match.arg(approximation)
  wilks <- mmemmuris:::lambda.mlm(fMod, within)

  if(approximation == "Rao"){
    raoFApprox <- mmemmuris:::raoF(wilks)

    withinTests <- raoFApprox$wilksWithin

    if(nrow(raoFApprox$wilksBetween$wilks) > 0L){
      betweenTests <- raoFApprox$wilksBetween
      wilks <- list(betweenTests = betweenTests,
                    withinTests = withinTests,
                    LR = FALSE)
      class(wilks) <- c("list", "multiTest", "wilks", "mlm")
    }else{
      wilks <- list(withinTests = withinTests,
                    LR = FALSE)
      class(wilks) <- c("list", "multiTest", "wilks", "mlm")
    }
    return(wilks)
  }else if(approximation == "LR"){
    lrChiApprox <- mmemmuris:::lrChi(wilks)

    withinTests <- lrChiApprox$wilksWithin

    if(nrow(lrChiApprox$wilksBetween$wilks) > 0L){
      betweenTests <- lrChiApprox$wilksBetween
      wilks <- list(betweenTests = betweenTests,
                    withinTests = withinTests,
                    LR = FALSE)
      class(wilks) <- c("list", "multiTest", "wilks", "mlm")
    }else{
      wilks <- list(withinTests = withinTests,
                    LR = FALSE)
      class(wilks) <- c("list", "multiTest", "wilks", "mlm")
    }
    return(wilks)
  }else if(approximation == "Bartlett"){
    BChiApprox <- mmemmuris:::BChi(wilks)

    withinTests <- BChiApprox$wilksWithin

    if(nrow(BChiApprox$wilksBetween$wilks) > 0L){
      betweenTests <- BChiApprox$wilksBetween
      wilks <- list(betweenTests = betweenTests,
                    withinTests = withinTests,
                    LR = FALSE)
      class(wilks) <- c("list", "multiTest", "wilks", "mlm")
    }else{
      wilks <- list(withinTests = withinTests,
                    LR = FALSE)
      class(wilks) <- c("list", "multiTest", "wilks", "mlm")
    }
    return(wilks)
  }
}

#' Multivariate Tests: Wilks' Lambda
#'
#' @export
#'
#'
#' @description
#' This function will calculate Wilks' \eqn{\Lambda} statistics from marginal models,
#' mixed-effects models, and multivariate linear models.
#'
#' @details
#'
#' **H and E matrices**
#'
#' The Wilks' lambda statistic for RM ANOVA is calculated as
#'
#' \eqn{\Lambda = (det(}**E**\eqn{)) / (det(}**H**\eqn{+}**E**\eqn{))}
#'
#' where the SSCP error matrix is
#'
#' **E = M'** \eqn{(}**Y'Y**\eqn{-} **\eqn{\hat{B}}'** \eqn{(}**X'X**\eqn{)}**\eqn{\hat{B}}** \eqn{)}**M**
#'
#' and the SSCP hypothesis matrix is
#'
#' **H** \eqn{=} **M'** \eqn{(}**L\eqn{\hat{B}}** \eqn{)}'\eqn{(}**L**\eqn{(}**X'X**\eqn{) ^ {-}}**L'** \eqn{) ^ {-1}}**L\eqn{\hat{B}}M**
#'
#' where **M** is an identity matrix for between-subject effects and a
#' sum-to-zero contrast matrix for within-subject effects, **Y** is a combined
#' matrix of observations, **\eqn{\hat{B}}** is a matrix of combined coefficients from the
#' ordinary least squares models, and **X** is the design matrix of the
#' independent variables.
#'
#' For marginal and mixed effects models, the **E** matrix is calculated from
#' the marginal variance-covariance matrix of the random/repeated effects
#' (**V**)
#'
#' **E =** \eqn{n}**M'VM**
#'
#' where **V = ZGZ' + R** (**Z** is the random effects design matrix, **G** is
#' the corresponding (co)variances of the random effects matrix, and **R**
#' represents the (co)variances of the repeated effects matrix).  The **H**
#' matrix is calculated with the same formula as above.  The **\eqn{\hat{B}}** matrix is
#' calculated by switching the reference levels of the within-subject effect,
#' refitting the models, and extracting the model coefficients.
#'
#' **Likelihood Ratio**
#'
#' The Wilks' lambda statistic for a marginal or mixed-effects model is
#' essentially a likelihood ratio test with the reduced model excluding the
#' term being tested.  The log-likelihood function for the full model is
#'
#' \eqn{logLik(M_F) = (-1 / 2)log(det(}**V_{F}** \eqn{)) - ((1 / 2)}**r'V_{F}** \eqn{^ {-1}}**r** \eqn{) - (N / 2)log(2\pi)}
#'
#' where **V_{F}** is the marginal variance-covariance matrix of the random/repeated
#' effects, **r** is the vector of model residuals, and \eqn{N} is the number of
#' observations of data in the long format.
#'
#' The contrast coefficients matrix (**L**) is calculated from the **LU**
#' decomposition of the crossproducts matrix (**X_{F}'X_{F}**).  To
#' project the design matrix of the reduced model (**X_{R}**) onto the
#' column space of the design matrix of the full model, the design matrix of the
#' full model needs to be multiplied by the orthogonal complement of **L'**.
#' It follows that the log-likelihood function for the reduced model is
#'
#' \eqn{logLik(M_R) = (-1 / 2)log(det(}**V_{R}** \eqn{)) - ((1 / 2)}**r'V_{R}** \eqn{^ {-1}}**r** \eqn{) - (N / 2)log(2\pi)}
#'
#' where **V_{R}** is the marginal variance-covariance matrix of the random/repeated
#' effects, **r** is the vector of model residuals, and \eqn{N} is the number of
#' observations of data in the long format.
#'
#' Wilks' \eqn{\Lambda} is denoted as
#' \deqn{\Lambda = exp{(-[-2logLik(M_R) - {-2logLik(M_F)}]) / n}}
#' where \eqn{n} is the number of subjects.
#'
#' The likelihood ratio (LR) statistic approximately follows a \eqn{\chi ^ 2}
#' distribution \deqn{LR = -2logLik(M_R) - {-2logLik(M_F)}} with
#' \eqn{df_{M_R} - df_{M_F}} degrees of freedom.
#'
#' **Rao**
#'
#' Rao's approximation approximately follows an F-distribution
#' \deqn{F_\Lambda = [(1 - \Lambda ^ {1 / t}) / (\Lambda ^ {1 / t})][(rt - 2u) / (pq)]}
#' where \eqn{p = rank(}**H**\eqn{+}**E**\eqn{)}, \eqn{q = rank(}**L**\eqn{(}**X'X**\eqn{)^{-1}}**L'** \eqn{)}, \eqn{v}
#' is the residual degrees of freedom, \eqn{r = v - (p - q + 1) / 2},
#' \eqn{u = (pq - 2) / 4}, and
#'
#' if \eqn{p ^ 2 + q ^ 2 - 5 > 0}, \deqn{t = \sqrt{(p ^ 2q ^ 2 - 4) / (p ^ 2 + q ^ 2 - 5)}}
#' else \deqn{t = 1}
#'
#' with \eqn{pq} numerator degrees of freedom and \eqn{rt - 2u} denominator
#' degrees of freedom.
#'
#' If \eqn{s = min(p, q) \le 2}, the F-statistic is exact.
#'
#' **Bartlett-Corrected Likelihood Ratio**
#'
#' The Bartlett-corrected LR statistic approximately follows a \eqn{\chi ^ 2}
#' distribution \deqn{LR_B = - (((p - q + 1) / 2) - v)ln(\Lambda)} with
#' \eqn{pq} degrees of freedom.
#'
#' @param fMod A model of class \code{\link[nlme]{gls}}, \code{\link[nlme]{lme}}, or \code{mlm}.
#' @param individual The subject whose marginal variance-covariance matrix
#' should be used, where **V = ZGZ' + R** (**Z** is the random effects design
#' matrix, **G** is the corresponding (co)variances of the random effects matrix,
#' and **R** represents the (co)variances of the repeated effects matrix).  The
#' default is the first subject in the factor.
#' @param approximation Approximation used for Wilks' \eqn{\Lambda}.  Options are Rao F ("Rao"),
#' Bartlett-corrected likelihood ratio \eqn{\chi ^ 2} ("Bartlett"), and likelihood
#' ratio \eqn{\chi ^ 2} ("LR").
#' @param within The name that should be used for the within-subject effect.
#' The default is "Time".
#' @param LR Likelihood ratio tests are used to calculate the likelihood ratio
#' \eqn{\chi ^ 2} and Wilks' \eqn{\Lambda} statistics instead of the **H** and
#' **E** matrices for marginal and mixed-effects models.
#'
#' @references {\url{https://support.sas.com/resources/papers/proceedings/proceedings/sugi23/Stats/p229.pdf}}
#' @references {\url{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.4/statug/statug_introreg_sect038.htm}}
#' @references {\url{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.3/statug/statug_mixed_details01.htm#statug.mixed.mixedmmtgandr}}
#' @references {Rao, C. R. (1973). Linear Statistical Inference and Its Applications. 2nd ed. New York: John Wiley & Sons.}
#'
#' @returns
#' This function will return two lists: a within-subject effects list
#' (`withinTests`) and a between-subject effects list (`betweenTests`) if
#' between-subject effects are in the model.
#' `withinTests` is a list containing: \tabular{ll}{
#'    `wilks`  \tab A `data.frame` containing the Wilks' \eqn{\Lambda} statistics and
#'    hypothesis tests for these statistics. \cr
#'    `parmsWilks`  \tab A `data.frame` containing the values used to calculate
#'    the Wilks' \eqn{\Lambda} statistics. \cr
#'    `approximation`  \tab Approximation used for the Wilks' \eqn{\Lambda} statistics.  Options are Rao F ("Rao"),
#' Bartlett-corrected likelihood ratio \eqn{\chi ^ 2} ("Bartlett"), and likelihood
#' ratio \eqn{\chi ^ 2} ("LR"). \cr
#'    `E`  \tab The sums of squares and crossproducts error matrix for within-subject effects. \cr
#'    `H`  \tab The sums of squares and crossproducts hypothesis matrices for within-subject effects. \cr
#'    `L`  \tab The contrast coefficients matrices. \cr
#'    `type`  \tab The type of **L** matrix used.  Type 1 is sequential. \cr
#' \tab \cr
#' }
#' `betweenTests` is a list containing: \tabular{ll}{
#'    `wilks`  \tab A `data.frame` containing the Wilks' \eqn{\Lambda} statistics and
#'    hypothesis tests for these statistics. \cr
#'    `parmsWilks`  \tab A `data.frame` containing the values used to calculate
#'    the Wilks' \eqn{\Lambda} statistics. \cr
#'    `approximation`  \tab Approximation used for the Wilks' \eqn{\Lambda} statistics.  Options are Rao F ("Rao"),
#' Bartlett-corrected likelihood ratio \eqn{\chi ^ 2} ("Bartlett"), and likelihood
#' ratio \eqn{\chi ^ 2} ("LR"). \cr
#'    `EB`  \tab The sums of squares and crossproducts error matrix for between-subject effects. \cr
#'    `HB`  \tab The sums of squares and crossproducts hypothesis matrices for between-subject effects. \cr
#'    `L`  \tab The contrast coefficients matrices. \cr
#'    `type`  \tab The type of **L** matrix used.  Type 1 is sequential. \cr
#'    `Vmatrix`  \tab The marginal variance-covariance matrix of the random/repeated effects. \cr
#' \tab \cr
#' }
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#' # Rao F statistic
#' wilksLambda(glsMod)
#' # Bartlett-corrected likelihood ratio chi-square statistic
#' wilksLambda(glsMod, approximation = "B")
#' # Likelihood ratio chi-square statistic
#' wilksLambda(glsMod, approximation = "LR")
#'
#' # Linear mixed-effects model
#' # Unstructured covariance
#' lmeMod <- lme(Distance ~ Sex * Age,
#'               random = ~ Age | Subject,
#'               na.action = na.omit,
#'               data = orthodontLong)
#' # Rao F statistic
#' wilksLambda(lmeMod)
#' # Bartlett-corrected likelihood ratio chi-square statistic
#' wilksLambda(lmeMod, approximation = "B")
#' # Likelihood ratio chi-square statistic
#' wilksLambda(lmeMod, approximation = "LR")
#'
#' # Multivariate linear model
#' manovaMod <- manova(cbind(Distance8, Distance10, Distance12, Distance14) ~ Sex,
#'                     data = orthodontWide)
#' # Rao F statistic
#' wilksLambda(manovaMod)
#' # Bartlett-corrected likelihood ratio chi-square statistic
#' wilksLambda(manovaMod, approximation = "B")
#' # Likelihood ratio chi-square statistic
#' wilksLambda(manovaMod, approximation = "LR")
#'
#' @seealso
#' {\code{\link[mmemmuris]{Hmatrix}}, \code{\link[mmemmuris]{coefs}}, \code{\link[mmemmuris]{Vmatrix}}, \code{\link[mmemmuris]{Ematrix}}, \code{\link[mmemmuris]{REMLtoML}}, \code{\link[mmemmuris]{LRT}}}

wilksLambda <- function(fMod, individual = NULL, approximation = c("Rao", "Bartlett", "LR"),
                        within = "Time", LR = TRUE, ss = NULL){
  if(any(class(fMod) %in% c("gls", "lme"))){
    wilks <- mmemmuris:::wilksLambda.ulm(fMod, individual, approximation, LR = LR, ss = ss)
    return(wilks)
  }else if(any(class(fMod) %in% "mlm")){
    wilks <- mmemmuris:::wilksLambda.mlm(fMod, approximation, within)
    return(wilks)
  }else
    stop('Please provide a model of class "gls", "lme", or "mlm".', call. = FALSE)
}

lambda.mlm <- function(fMod, within = "Time"){
  type <- "1"
  E <- mmemmuris:::Ematrix.mlm(fMod)
  n <- E$n
  H <- mmemmuris:::Hmatrix.mlm(fMod, within)

  # Calculate Wilks' lambda statistic
  lambdaBetween <- sapply(H$HB, function(x){
    det(E$EB) / det(x + E$EB)
  })
  lambdaWithin <- sapply(H$H, function(x){
    det(E$E) / det(x + E$E)
  })

  pBetween <- sapply(H$HB, function(x){
    pBetween <- Matrix::rankMatrix(x + E$EB)
    attributes(pBetween) <- NULL
    pBetween
  })
  pWithin <- sapply(H$H, function(x){
    pWithin <- Matrix::rankMatrix(x + E$E)
    attributes(pWithin) <- NULL
    pWithin
  })

  X <- H$X
  Qbetween <- Qwithin <- sapply(H$L, function(x){
    Q <- Matrix::rankMatrix(x %*% MASS::ginv(t(X) %*% X) %*% t(x))
    attributes(Q) <- NULL
    Q
  })

  rankX <- Matrix::rankMatrix(X)
  attributes(rankX) <- NULL
  v <- n - rankX

  sBetween <- apply(cbind(pBetween, Qbetween), 1, min)
  sWithin <- apply(cbind(pWithin, Qwithin), 1, min)

  X2Between <- -n * log(lambdaBetween)
  X2Within <- -n * log(lambdaWithin)

  lambda <- list(EB = E$EB, HB = H$HB, lambdaBetween = lambdaBetween,
                 sBetween = sBetween, pBetween = pBetween, qBetween = Qbetween,
                 X2Between = X2Between, E = E$E, H = H$H,
                 lambdaWithin = lambdaWithin, sWithin = sWithin,
                 pWithin = pWithin, qWithin = Qwithin, X2Within = X2Within, X = X, v = v,
                 Lbetween = H$L, n = n, type = type)

  class(lambda) <- c("list", "lambda.mlm")
  return(lambda)
}

lambda.ulm <- function(fMod, individual = NULL, LR = TRUE, ss = NULL){
  type <- "1"
  if(mmemmuris::covStruct(fMod) %in% c("cs", "other"))
    stop("Only unstructured covariance models are allowed.", call. = FALSE)
  tt <- mmemmuris::termsType(fMod)
  HBsingular <- Hsingular <- FALSE
  E <- mmemmuris:::Ematrix.ulm(fMod, individual, ss)
  H <- mmemmuris:::Hmatrix.ulm(fMod)
  n <- E$n
  if("error" %in% class(tryCatch(any(sapply(H$HB, function(x){ mmemmuris:::determinantMatrix(x + E$EB) == 0L })), error = function(e) { e })))
    stop('Please pick a different E matrix with the "individual" argument.', call. = FALSE)
  if("error" %in% class(tryCatch(any(sapply(H$H, function(x){ mmemmuris:::determinantMatrix(x + E$E) == 0L })), error = function(e) { e })))
    stop('Please pick a different E matrix with the "individual" argument.', call. = FALSE)
  if(any(sapply(H$HB, function(x){ mmemmuris:::determinantMatrix(x + E$EB) == 0L })))
    HBsingular <- TRUE
  if(any(sapply(H$H, function(x){ mmemmuris:::determinantMatrix(x + E$E) == 0L })))
    Hsingular <- TRUE
  if(HBsingular == TRUE | Hsingular == TRUE){
    LR <- TRUE
    warning("det(HB + EB) = 0 or det(H + E) = 0.  Likelihood ratio tests used instead.", call. = FALSE)
  }
  if(LR == FALSE){
    lambdaBetween <- sapply(H$HB, function(x){
      det(E$EB) / det(x + E$EB)
    })
    lambdaWithin <- sapply(H$H, function(x){
      det(E$E) / det(x + E$E)
    })

    pBetween <- sapply(H$HB, function(x){
      pBetween <- Matrix::rankMatrix(x + E$EB)
      attributes(pBetween) <- NULL
      pBetween
    })
    pWithin <- sapply(H$H, function(x){
      pWithin <- Matrix::rankMatrix(x + E$E)
      attributes(pWithin) <- NULL
      pWithin
    })

    X <- H$X
    Lbetween <- H$Lbetween
    L <- Lwithin <- NULL
  }else{
    HB <- H <- NULL
    L <- mmemmuris:::Lmatrix.ulm(fMod)
    Lwithin <- L[tt$within]
    if(length(tt$between) == 0L){
      L <- LRTbetween <- X2Between <- lambdaBetween <- NULL
    }else{
      term <- attr(terms(fMod), "term.labels")
      index <- grepl(paste0(":", tt$within[1], "|", tt$within[1], ":"), term)
      index <- unname(as.list(data.frame(rbind(tt$between, term[index]))))
      L <- lapply(index, function(x){ do.call("rbind", mmemmuris:::Lmatrix.ulm(fMod)[x]) })
      names(L) <- tt$between
      LRTbetween <- mmemmuris::LRT(fMod, L)

      X2Between <- sapply(LRTbetween, function(x){
        x$L.Ratio[2]
      })

      lambdaBetween <- exp(-(X2Between) / n)
    }

    LRTwithin <- mmemmuris::LRT(fMod, Lwithin)

    X2Within <- sapply(LRTwithin, function(x){
      x$L.Ratio[2]
    })

    lambdaWithin <- exp(-(X2Within) / n)

    X <- tt$Xbetween

    Lbetween <- t(as.matrix(Matrix::expand(Matrix::lu(t(X) %*% X))$L))
    rownames(Lbetween) <- colnames(Lbetween) <- colnames(X)

    termsIndex <- attr(X, "assign")
    tableTerms <- table(termsIndex)
    termsLabels <- rep(c("(Intercept)", attr(terms(fMod), "term.labels")[attr(terms(fMod), "term.labels") %in% tt$between]), tableTerms)
    term <- c("(Intercept)", attr(terms(fMod), "term.labels")[attr(terms(fMod), "term.labels") %in% tt$between])
    rowIndex <- list()
    for(i in 1:length(term)){
      rowIndex[[i]] <- which(termsLabels == term[i])
      names(rowIndex)[i] <- term[i]
    }
    Lbetween <- lapply(rowIndex, function(x){
      Lbetween[x,, drop = FALSE]
    })
  }
  Qbetween <- Qwithin <- sapply(Lbetween, function(x){
    Q <- Matrix::rankMatrix(x %*% MASS::ginv(t(X) %*% X) %*% t(x))
    attributes(Q) <- NULL
    Q
  })

  rankX <- Matrix::rankMatrix(X)
  attributes(rankX) <- NULL
  v <- n - rankX

  if(LR == FALSE){
    sBetween <- apply(cbind(pBetween, Qbetween), 1, min)
    X2Between <- -n * log(lambdaBetween)
    X2Within <- -n * log(lambdaWithin)
    Qwithin <- Qwithin[1:length(tt$within)]
  }else{
    if(length(tt$between) == 0L){
      Qbetween <- pBetween <- sBetween <- NULL
    }else{
      pqBetween <- sapply(LRTbetween, function(x){
        x$df[2] - x$df[1]
      })
      Qbetween <- Qbetween[-1]
      pBetween <- pqBetween / Qbetween
      #names(pBetween) <- tt$between
      sBetween <- apply(cbind(pBetween, Qbetween), 1, min)
    }
    pqWithin <- sapply(LRTwithin, function(x){
      x$df[2] - x$df[1]
    })
    Qwithin <- Qwithin[1:length(tt$within)]
    pWithin <- pqWithin / Qwithin
    names(pWithin) <- tt$within
  }
  sWithin <- apply(cbind(pWithin, Qwithin), 1, min)[1:length(tt$within)]

  lambda <- list(EB = E$EB, HB = H$HB, lambdaBetween = lambdaBetween,
                 sBetween = sBetween, pBetween = pBetween, qBetween = Qbetween,
                 X2Between = X2Between, E = E$E, H = H$H,
                 lambdaWithin = lambdaWithin, sWithin = sWithin,
                 pWithin = pWithin, qWithin = Qwithin, X2Within = X2Within,
                 X = X, v = v, Lbetween = Lbetween, Lwithin = Lwithin, L = L, n = n, type = type, LR = LR)

  class(lambda) <- c("list", "lambda.ulm")
  return(lambda)
}

V.mlm <- function(fMod, within = "Time"){
  type <- "1"
  E <- mmemmuris:::Ematrix.mlm(fMod)
  n <- E$n
  H <- mmemmuris:::Hmatrix.mlm(fMod, within)

  # Calculate Pillai's trace statistic
  VBetween <- sapply(H$HB, function(x){
    sum(diag(x %*% solve(x + E$EB)))
  })
  VWithin <- sapply(H$H, function(x){
    sum(diag(x %*% solve(x + E$E)))
  })

  pBetween <- sapply(H$HB, function(x){
    pBetween <- Matrix::rankMatrix(x + E$EB)
    attributes(pBetween) <- NULL
    pBetween
  })
  pWithin <- sapply(H$H, function(x){
    pWithin <- Matrix::rankMatrix(x + E$E)
    attributes(pWithin) <- NULL
    pWithin
  })

  X <- H$X
  Qbetween <- Qwithin <- sapply(H$L, function(x){
    Q <- Matrix::rankMatrix(x %*% MASS::ginv(t(X) %*% X) %*% t(x))
    attributes(Q) <- NULL
    Q
  })

  rankX <- Matrix::rankMatrix(X)
  attributes(rankX) <- NULL
  v <- n - rankX

  sBetween <- apply(cbind(pBetween, Qbetween), 1, min)
  sWithin <- apply(cbind(pWithin, Qwithin), 1, min)

  V <- list(EB = E$EB, HB = H$HB, VBetween = VBetween,
            sBetween = sBetween, pBetween = pBetween, qBetween = Qbetween,
            E = E$E, H = H$H,
            VWithin = VWithin, sWithin = sWithin,
            pWithin = pWithin, qWithin = Qwithin, X = X, v = v,
            Lbetween = H$L, n = n, type = type)

  class(V) <- c("list", "V.mlm")
  return(V)
}

V.ulm <- function(fMod, individual = NULL, LR = TRUE, ss = NULL){
  type <- "1"
  if(mmemmuris::covStruct(fMod) %in% c("cs", "other"))
    stop("Only unstructured covariance models are allowed.", call. = FALSE)
  tt <- mmemmuris::termsType(fMod)
  HBsingular <- Hsingular <- FALSE
  E <- mmemmuris:::Ematrix.ulm(fMod, individual, ss)
  H <- mmemmuris:::Hmatrix.ulm(fMod)
  n <- E$n
  if("error" %in% class(tryCatch(any(sapply(H$HB, function(x){ is.null(mmemmuris:::inverseMatrix(x + E$EB)) })), error = function(e) { e })))
    stop('Please pick a different E matrix with the "individual" argument.', call. = FALSE)
  if("error" %in% class(tryCatch(any(sapply(H$H, function(x){ is.null(mmemmuris:::inverseMatrix(x + E$E)) })), error = function(e) { e })))
    stop('Please pick a different E matrix with the "individual" argument.', call. = FALSE)
  if(any(sapply(H$HB, function(x){ is.null(mmemmuris:::inverseMatrix(x + E$EB)) })))
    HBsingular <- TRUE
  if(any(sapply(H$H, function(x){ is.null(mmemmuris:::inverseMatrix(x + E$E)) })))
    Hsingular <- TRUE
  if(HBsingular == TRUE | Hsingular == TRUE){
    LR <- TRUE
    warning("(HB + EB) ^ {-1} or (H + E) ^ {-1} are singular matrices.  Likelihood ratio tests (s = 1) used instead.", call. = FALSE)
  }
  if(LR == FALSE){
    VBetween <- sapply(H$HB, function(x){
      sum(diag(x %*% solve(x + E$EB)))
    })
    VWithin <- sapply(H$H, function(x){
      sum(diag(x %*% solve(x + E$E)))
    })

    pBetween <- sapply(H$HB, function(x){
      pBetween <- Matrix::rankMatrix(x + E$EB)
      attributes(pBetween) <- NULL
      pBetween
    })
    pWithin <- sapply(H$H, function(x){
      pWithin <- Matrix::rankMatrix(x + E$E)
      attributes(pWithin) <- NULL
      pWithin
    })

    X <- H$X
    Lbetween <- H$Lbetween
    L <- Lwithin <- NULL
  }else{
    HB <- H <- NULL
    L <- mmemmuris:::Lmatrix.ulm(fMod)
    Lwithin <- L[tt$within]
    if(length(tt$between) == 0L){
      L <- LRTbetween <- X2Between <- VBetween <- NULL
    }else{
      term <- attr(terms(fMod), "term.labels")
      index <- grepl(paste0(":", tt$within[1], "|", tt$within[1], ":"), term)
      index <- unname(as.list(data.frame(rbind(tt$between, term[index]))))
      L <- lapply(index, function(x){ do.call("rbind", mmemmuris:::Lmatrix.ulm(fMod)[x]) })
      names(L) <- tt$between
      LRTbetween <- mmemmuris::LRT(fMod, L)

      X2Between <- sapply(LRTbetween, function(x){
        x$L.Ratio[2]
      })

      VBetween <- 1 - exp(-(X2Between) / n)
    }

    LRTwithin <- mmemmuris::LRT(fMod, Lwithin)

    X2Within <- sapply(LRTwithin, function(x){
      x$L.Ratio[2]
    })

    VWithin <- 1 - exp(-(X2Within) / n)

    X <- tt$Xbetween

    Lbetween <- t(as.matrix(Matrix::expand(Matrix::lu(t(X) %*% X))$L))
    rownames(Lbetween) <- colnames(Lbetween) <- colnames(X)

    termsIndex <- attr(X, "assign")
    tableTerms <- table(termsIndex)
    termsLabels <- rep(c("(Intercept)", attr(terms(fMod), "term.labels")[attr(terms(fMod), "term.labels") %in% tt$between]), tableTerms)
    term <- c("(Intercept)", attr(terms(fMod), "term.labels")[attr(terms(fMod), "term.labels") %in% tt$between])
    rowIndex <- list()
    for(i in 1:length(term)){
      rowIndex[[i]] <- which(termsLabels == term[i])
      names(rowIndex)[i] <- term[i]
    }
    Lbetween <- lapply(rowIndex, function(x){
      Lbetween[x,, drop = FALSE]
    })
  }
  Qbetween <- Qwithin <- sapply(Lbetween, function(x){
    Q <- Matrix::rankMatrix(x %*% MASS::ginv(t(X) %*% X) %*% t(x))
    attributes(Q) <- NULL
    Q
  })

  rankX <- Matrix::rankMatrix(X)
  attributes(rankX) <- NULL
  v <- n - rankX

  if(LR == FALSE){
    sBetween <- apply(cbind(pBetween, Qbetween), 1, min)
    Qwithin <- Qwithin[1:length(tt$within)]
  }else{
    if(length(tt$between) == 0L){
      Qbetween <- pBetween <- sBetween <- NULL
    }else{
      pqBetween <- sapply(LRTbetween, function(x){
        x$df[2] - x$df[1]
      })
      Qbetween <- Qbetween[-1]
      pBetween <- pqBetween / Qbetween
      #names(pBetween) <- tt$between
      sBetween <- apply(cbind(pBetween, Qbetween), 1, min)
    }
    pqWithin <- sapply(LRTwithin, function(x){
      x$df[2] - x$df[1]
    })
    Qwithin <- Qwithin[1:length(tt$within)]
    pWithin <- pqWithin / Qwithin
    names(pWithin) <- tt$within
  }
  sWithin <- apply(cbind(pWithin, Qwithin), 1, min)[1:length(tt$within)]

  V <- list(EB = E$EB, HB = H$HB, VBetween = VBetween,
            sBetween = sBetween, pBetween = pBetween, qBetween = Qbetween,
            E = E$E, H = H$H,
            VWithin = VWithin, sWithin = sWithin,
            pWithin = pWithin, qWithin = Qwithin, X = X, v = v,
            Lbetween = Lbetween, Lwithin = Lwithin, L = L, n = n, type = type, LR = LR)

  class(V) <- c("list", "V.ulm")
  return(V)
}

pillaiTrace.mlm <- function(fMod, approximation = c("Muller", "Pillai"), within = "Time"){

  approximation <- match.arg(approximation)
  V <- mmemmuris:::V.mlm(fMod, within)

  if(approximation == "Muller"){
    mullerFApprox <- mmemmuris:::mullerF(V)

    withinTests <- mullerFApprox$pillaiWithin

    if(nrow(mullerFApprox$pillaiBetween$pillai) > 0L){
      betweenTests <- mullerFApprox$pillaiBetween
      pillai <- list(betweenTests = betweenTests,
                     withinTests = withinTests,
                     LR = FALSE)
      class(pillai) <- c("list", "multiTest", "pillai", "mlm")
    }else{
      pillai <- list(withinTests = withinTests,
                     LR = FALSE)
      class(pillai) <- c("list", "multiTest", "pillai", "mlm")
    }
    return(pillai)
  }else{
    pillaiFApprox <- mmemmuris:::pillaiF(V)

    withinTests <- pillaiFApprox$pillaiWithin

    if(nrow(pillaiFApprox$pillaiBetween$pillai) > 0L){
      betweenTests <- pillaiFApprox$pillaiBetween
      pillai <- list(betweenTests = betweenTests,
                     withinTests = withinTests,
                     LR = FALSE)
      class(pillai) <- c("list", "multiTest", "pillai", "mlm")
    }else{
      pillai <- list(withinTests = withinTests,
                     LR = FALSE)
      class(pillai) <- c("list", "multiTest", "pillai", "mlm")
    }
    return(pillai)
  }
}

pillaiTrace.ulm <- function(fMod, individual = NULL, approximation = c("Muller", "Pillai"), LR = TRUE, ss = NULL){
  approximation <- match.arg(approximation)
  tt <- mmemmuris::termsType(fMod)
  V <- mmemmuris:::V.ulm(fMod, individual, LR = LR, ss = ss)

  if(approximation == "Muller"){
    mullerFApprox <- mmemmuris:::mullerF(V)

    withinTests <- mullerFApprox$pillaiWithin

    if(!is.null(mullerFApprox$pillaiBetween$pillai)){
      if(nrow(mullerFApprox$pillaiBetween$pillai) > 0L){
        betweenTests <- mullerFApprox$pillaiBetween
        pillai <- list(betweenTests = betweenTests,
                       withinTests = withinTests,
                       LR = V$LR)
        class(pillai) <- c("list", "multiTest", "pillai", "ulm")
      }else{
        pillai <- list(withinTests = withinTests,
                       LR = V$LR)
        class(pillai) <- c("list", "multiTest", "pillai", "ulm")
      }
    }else{
      pillai <- list(withinTests = withinTests,
                     LR = V$LR)
      class(pillai) <- c("list", "multiTest", "pillai", "ulm")
    }
    return(pillai)
  }else if(approximation == "Pillai"){
    pillaiFApprox <- mmemmuris:::pillaiF(V)

    withinTests <- pillaiFApprox$pillaiWithin

    if(!is.null(pillaiFApprox$pillaiBetween$pillai)){
      if(nrow(pillaiFApprox$pillaiBetween$pillai) > 0L){
        betweenTests <- pillaiFApprox$pillaiBetween
        pillai <- list(betweenTests = betweenTests,
                       withinTests = withinTests,
                       LR = V$LR)
        class(pillai) <- c("list", "multiTest", "pillai", "ulm")
      }else{
        pillai <- list(withinTests = withinTests,
                       LR = V$LR)
        class(pillai) <- c("list", "multiTest", "pillai", "ulm")
      }
    }else{
      pillai <- list(withinTests = withinTests,
                     LR = V$LR)
      class(pillai) <- c("list", "multiTest", "pillai", "ulm")
    }
    return(pillai)
  }
}

#' Multivariate Tests: Pillai's Trace
#'
#' @export
#'
#' @description
#' This function will calculate Pillai's Trace \eqn{V} statistics from marginal models,
#' mixed-effects models, and multivariate linear models.
#'
#' @details
#'
#' **H and E matrices**
#'
#' The Pillai's trace statistic for a RM ANOVA is denoted as
#'
#' \eqn{V = trace(}**H**\eqn{(}**H**\eqn{+}**E**\eqn{)^{-1})}
#'
#' where the SSCP error matrix is
#'
#' **E = M'** \eqn{(}**Y'Y**\eqn{-} **\eqn{\hat{B}}'** \eqn{(}**X'X**\eqn{)}**\eqn{\hat{B}}** \eqn{)}**M**
#'
#' and the SSCP hypothesis matrix is
#'
#' **H** \eqn{=} **M'** \eqn{(}**L\eqn{\hat{B}}** \eqn{)}'\eqn{(}**L**\eqn{(}**X'X**\eqn{) ^ {-}}**L'** \eqn{) ^ {-1}}**L\eqn{\hat{B}}M**
#'
#' where **M** is an identity matrix for between-subject effects and a
#' sum-to-zero contrast matrix for within-subject effects, **Y** is a combined
#' matrix of observations, **\eqn{\hat{B}}** is a matrix of combined coefficients from the
#' ordinary least squares models, and **X** is the design matrix of the
#' independent variables.
#'
#' For marginal and mixed effects models, the **E** matrix is calculated from
#' the marginal variance-covariance matrix of the random/repeated effects
#' (**V**)
#'
#' **E =** \eqn{n}**M'VM**
#'
#' where **V = ZGZ' + R** (**Z** is the random effects design matrix, **G** is
#' the corresponding (co)variances of the random effects matrix, and **R**
#' represents the (co)variances of the repeated effects matrix).  The **H**
#' matrix is calculated with the same formula as above.  The **\eqn{\hat{B}}** matrix is
#' calculated by switching the reference levels of the within-subject effect,
#' refitting the models, and extracting the model coefficients.
#'
#' **Likelihood Ratio**
#'
#' The Wilks' lambda statistic for a marginal or mixed-effects model is
#' essentially a likelihood ratio test with the reduced model excluding the
#' term being tested.  The log-likelihood function for the full model is
#'
#' \eqn{logLik(M_F) = (-1 / 2)log(det(}**V_{F}** \eqn{)) - ((1 / 2)}**r'V_{F}** \eqn{^ {-1}}**r** \eqn{) - (N / 2)log(2\pi)}
#'
#' where **V_{F}** is the marginal variance-covariance matrix of the random/repeated
#' effects, **r** is the vector of model residuals, and \eqn{N} is the number of
#' observations of data in the long format.
#'
#' The contrast coefficients matrix (**L**) is calculated from the **LU**
#' decomposition of the crossproducts matrix (**X_{F}'X_{F}**).  To
#' project the design matrix of the reduced model (**X_{R}**) onto the
#' column space of the design matrix of the full model, the design matrix of the
#' full model needs to be multiplied by the orthogonal complement of **L'**.
#' It follows that the log-likelihood function for the reduced model is
#'
#' \eqn{logLik(M_R) = (-1 / 2)log(det(}**V_{R}** \eqn{)) - ((1 / 2)}**r'V_{R}** \eqn{^ {-1}}**r** \eqn{) - (N / 2)log(2\pi)}
#'
#' where **V_{R}** is the marginal variance-covariance matrix of the random/repeated
#' effects, **r** is the vector of model residuals, and \eqn{N} is the number of
#' observations of data in the long format.
#'
#' Wilks' \eqn{\Lambda} is denoted as
#' \deqn{\Lambda = exp{(-[-2logLik(M_R) - {-2logLik(M_F)}]) / n}}
#' where \eqn{n} is the number of subjects.
#'
#' For marginal and mixed-effects models, if \eqn{s = min(p, q) = 1}, Pillai's
#' trace is calculated as \eqn{V = 1 - \Lambda}.  With `LR = TRUE`,
#' Pillai's trace is only calculated when \eqn{s = 1} for within-subject
#' effects.
#'
#' **Pillai (1954)**
#'
#' Let \eqn{v} be the residual degrees of freedom, \eqn{p = rank(}**H**\eqn{+}**E**\eqn{)},
#' \eqn{q = rank(}**L**\eqn{(}**X'X**\eqn{)^{-1}}**L'** \eqn{)}, \eqn{s = min(p, q)},
#' \eqn{m = (|p - q| - 1) / 2} and \eqn{n = (v - p - 1) / 2}.  The
#' statistic approximately follows an F-distribution
#' \deqn{F_V = ((2n + s + 1) / (2m + s + 1))(V / (s - V))} with
#' \eqn{s(2m + s + 1)} numerator degrees of freedom and \eqn{s(2n + s + 1)}
#' denominator degrees of freedom.
#'
#' If \eqn{s = 1}, the F-statistic is exact.
#'
#' **Muller Method #2 (1998)**
#'
#' Let \eqn{v} be the residual degrees of freedom, \eqn{p = rank(}**H**\eqn{+}**E**\eqn{)},
#' \eqn{q = rank(}**L**\eqn{(}**X'X**\eqn{)^{-1}}**L'** \eqn{)}, and \eqn{s = min(p, q)}.  The
#' statistic approximately follows an F-distribution
#' \deqn{F_V = ((((v + s - p) / (v + q))[(((s(v + s - p)(v + q + 2)(v + q - 1)) / (v(v + q - p))) - 2)]) / (((qp) / (s(v + q)))[(((s(v + s - p)(v + q + 2)(v + q - 1)) / (v(v + q - p))) - 2)]))(V / (s - V))} with
#' \eqn{((qp) / (s(v + q)))[(((s(v + s - p)(v + q + 2)(v + q - 1)) / (v(v + q - p))) - 2)]} numerator degrees of freedom and \eqn{(((v + s - p) / (v + q))[(((s(v + s - p)(v + q + 2)(v + q - 1)) / (v(v + q - p))) - 2)])}
#' denominator degrees of freedom.
#'
#' If \eqn{s = 1}, the F-statistic is exact.
#'
#' @param fMod A model of class \code{\link[nlme]{gls}}, \code{\link[nlme]{lme}}, or \code{mlm}.
#' @param individual The subject whose marginal variance-covariance matrix
#' should be used, where **V = ZGZ' + R** (**Z** is the random effects design
#' matrix, **G** is the corresponding (co)variances of the random effects matrix,
#' and **R** represents the (co)variances of the repeated effects matrix).  The
#' default is the first subject in the factor.
#' @param within The name that should be used for the within-subject effect.
#' The default is "Time".
#' @param approximation Approximation used for Pillai's trace.  Options are
#' Pillai F ("Pillai") and Muller (Method #2) F ("Muller").
#' @param LR Likelihood ratio tests are used to calculate the likelihood ratio
#' \eqn{\chi ^ 2} and Wilks' \eqn{\Lambda} statistics instead of the **H** and
#' **E** matrices.
#'
#' @returns
#' This function will return two lists: a within-subject effects list
#' (`withinTests`) and a between-subject effects list (`betweenTests`) if
#' between-subject effects are in the model.
#' @returns `withinTests` is a list containing: \tabular{ll}{
#'    `pillai`  \tab A `data.frame` containing the Pillai's Trace \eqn{V} statistics and
#'    hypothesis tests for these statistics. \cr
#'    `parmsPillai`  \tab A `data.frame` containing the values used to calculate
#'    the Pillai's Trace \eqn{V} statistics. \cr
#'    `approximation`  \tab Approximation used for the Pillai's Trace \eqn{V} statistics.  Options are Muller F ("Muller") and
#' and Pillai F ("Pillai"). \cr
#'    `E`  \tab The sums of squares and crossproducts error matrix for within-subject effects. \cr
#'    `H`  \tab The sums of squares and crossproducts hypothesis matrices for within-subject effects. \cr
#'    `L`  \tab The contrast coefficients matrices. \cr
#'    `type`  \tab The type of **L** matrix used.  Type 1 is sequential. \cr
#' \tab \cr
#' }
#' `betweenTests` is a list containing: \tabular{ll}{
#'    `pillai`  \tab A `data.frame` containing the Pillai's Trace \eqn{V} statistics and
#'    hypothesis tests for these statistics. \cr
#'    `parmsPillai`  \tab A `data.frame` containing the values used to calculate
#'    the Pillai's Trace \eqn{V} statistics. \cr
#'    `approximation`  \tab Approximation used for the Pillai's Trace \eqn{V} statistics.  Options are Muller F ("Muller") and
#' and Pillai F ("Pillai"). \cr
#'    `EB`  \tab The sums of squares and crossproducts error matrix for between-subject effects. \cr
#'    `HB`  \tab The sums of squares and crossproducts hypothesis matrices for between-subject effects. \cr
#'    `L`  \tab The contrast coefficients matrices. \cr
#'    `type`  \tab The type of **L** matrix used.  Type 1 is sequential. \cr
#'    `Vmatrix`  \tab The marginal variance-covariance matrix of the random/repeated effects. \cr
#' \tab \cr
#' }
#'
#' @references {\url{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.4/statug/statug_introreg_sect038.htm}}
#' @references {Muller, K. E. (1998). "A New F Approximation for the Pillai-Bartlett Trace under \eqn{H_0}." Journal of Computational and Graphical Statistics 7:131-137.}
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#' # Pillai F statistic
#' pillaiTrace(glsMod)
#'
#' # Linear mixed-effects model
#' # Unstructured covariance
#' lmeMod <- lme(Distance ~ Sex * Age,
#'               random = ~ Age | Subject,
#'               na.action = na.omit,
#'               data = orthodontLong)
#' # Pillai F statistic
#' pillaiTrace(lmeMod)
#'
#' # Multivariate linear model
#' manovaMod <- manova(cbind(Distance8, Distance10, Distance12, Distance14) ~ Sex,
#'                     data = orthodontWide)
#' # Pillai F statistic
#' pillaiTrace(manovaMod)
#' # Muller F statistic (Method #2)
#' pillaiTrace(manovaMod, approximation = "M")
#'
#' @seealso
#' {\code{\link[mmemmuris]{Hmatrix}}, \code{\link[mmemmuris]{coefs}}, \code{\link[mmemmuris]{Vmatrix}}, \code{\link[mmemmuris]{Ematrix}}, \code{\link[mmemmuris]{REMLtoML}}, \code{\link[mmemmuris]{LRT}}}

pillaiTrace <- function(fMod, individual = NULL,
                        approximation = c("Muller", "Pillai"), within = "Time", LR = TRUE, ss = NULL){
  approximation <- match.arg(approximation)
  if(any(class(fMod) %in% c("gls", "lme"))){
    pillai <- mmemmuris:::pillaiTrace.ulm(fMod, individual, approximation, LR = LR, ss = ss)
    return(pillai)
  }else if(any(class(fMod) %in% "mlm")){
    pillai <- mmemmuris:::pillaiTrace.mlm(fMod, approximation, within)
    return(pillai)
  }else
    stop('Please provide a model of class "gls", "lme", or "mlm".', call. = FALSE)
}

U.mlm <- function(fMod, within = "Time"){
  type <- "1"
  E <- mmemmuris:::Ematrix.mlm(fMod)
  n <- E$n
  H <- mmemmuris:::Hmatrix.mlm(fMod, within)

  # Calculate Hotelling-Lawley trace statistic
  UBetween <- sapply(H$HB, function(x){
    sum(diag(solve(E$EB) %*% x))
  })
  UWithin <- sapply(H$H, function(x){
    sum(diag(solve(E$E) %*% x))
  })

  pBetween <- sapply(H$HB, function(x){
    pBetween <- Matrix::rankMatrix(x + E$EB)
    attributes(pBetween) <- NULL
    pBetween
  })
  pWithin <- sapply(H$H, function(x){
    pWithin <- Matrix::rankMatrix(x + E$E)
    attributes(pWithin) <- NULL
    pWithin
  })

  X <- H$X
  Qbetween <- Qwithin <- sapply(H$L, function(x){
    Q <- Matrix::rankMatrix(x %*% MASS::ginv(t(X) %*% X) %*% t(x))
    attributes(Q) <- NULL
    Q
  })

  rankX <- Matrix::rankMatrix(X)
  attributes(rankX) <- NULL
  v <- n - rankX

  X2Between <- v * UBetween
  X2Within <- v * UWithin

  sBetween <- apply(cbind(pBetween, Qbetween), 1, min)
  sWithin <- apply(cbind(pWithin, Qwithin), 1, min)

  U <- list(EB = E$EB, HB = H$HB, UBetween = UBetween,
            sBetween = sBetween, pBetween = pBetween, qBetween = Qbetween,
            X2Between = X2Between, E = E$E, H = H$H,
            UWithin = UWithin, sWithin = sWithin,
            pWithin = pWithin, qWithin = Qwithin, X2Within = X2Within, X = X, v = v,
            Lbetween = H$L, n = n, type = type)

  class(U) <- c("list", "U.mlm")
  return(U)
}

U.ulm <- function(fMod, individual = NULL, waldF = TRUE, ss = NULL){
  type <- "1"
  if(mmemmuris::covStruct(fMod) %in% c("cs", "other"))
    stop("Only unstructured covariance models are allowed.", call. = FALSE)
  tt <- mmemmuris::termsType(fMod)
  EBsingular <- Esingular <- FALSE
  E <- mmemmuris:::Ematrix.ulm(fMod, individual, ss)
  H <- mmemmuris:::Hmatrix.ulm(fMod)
  n <- E$n
  if("error" %in% class(tryCatch(any(sapply(H$HB, function(x){ Matrix::rankMatrix(x + E$EB) })), error = function(e) { e })))
    stop('Please pick a different E matrix with the "individual" argument.', call. = FALSE)
  if("error" %in% class(tryCatch(any(sapply(H$H, function(x){ Matrix::rankMatrix(x + E$E) })), error = function(e) { e })))
    stop('Please pick a different E matrix with the "individual" argument.', call. = FALSE)
  if(is.null(mmemmuris:::inverseMatrix(E$EB)))
    EBsingular <- TRUE
  if(is.null(mmemmuris:::inverseMatrix(E$E)))
    Esingular <- TRUE
  if(EBsingular == TRUE | Esingular == TRUE){
    waldF <- TRUE
    warning("EB ^ {-1} or E ^ {-1} are singular matrices.  Wald-type F-tests used instead.", call. = FALSE)
  }
  if(waldF == FALSE){
    UBetween <- sapply(H$HB, function(x){
      sum(diag(solve(E$EB) %*% x))
    })
    UWithin <- sapply(H$H, function(x){
      sum(diag(solve(E$E) %*% x))
    })

    pBetween <- sapply(H$HB, function(x){
      pBetween <- Matrix::rankMatrix(x + E$EB)
      attributes(pBetween) <- NULL
      pBetween
    })
    pWithin <- sapply(H$H, function(x){
      pWithin <- Matrix::rankMatrix(x + E$E)
      attributes(pWithin) <- NULL
      pWithin
    })

    X <- H$X
    rankX <- Matrix::rankMatrix(X)
    attributes(rankX) <- NULL
    v <- n - rankX
    Lbetween <- H$Lbetween
    L <- Lwithin <- NULL
  }else{
    HB <- H <- NULL
    X <- tt$Xbetween
    rankX <- Matrix::rankMatrix(X)
    attributes(rankX) <- NULL
    v <- n - rankX
    fMod <- mmemmuris::MLtoREML(fMod)

    L <- mmemmuris:::Lmatrix.ulm(fMod)
    Lwithin <- L[tt$within]
    if(length(tt$between) == 0L){
      L <- waldBetween <- X2Between <- UBetween <- NULL
    }else{
      term <- attr(terms(fMod), "term.labels")
      index <- grepl(paste0(":", tt$within[1], "|", tt$within[1], ":"), term)
      index <- unname(as.list(data.frame(rbind(tt$between, term[index]))))
      L <- lapply(index, function(x){ do.call("rbind", mmemmuris:::Lmatrix.ulm(fMod)[x]) })
      names(L) <- tt$between
      waldBetween <- mmemmuris::waldF(fMod, L)

      X2Between <- waldBetween$X2
      UBetween <- X2Between / v
    }

    waldWithin <- mmemmuris::waldF(fMod, Lwithin)
    X2Within <- waldWithin$X2
    UWithin <- X2Within / v

    Lbetween <- t(as.matrix(Matrix::expand(Matrix::lu(t(X) %*% X))$L))
    rownames(Lbetween) <- colnames(Lbetween) <- colnames(X)

    termsIndex <- attr(X, "assign")
    tableTerms <- table(termsIndex)
    termsLabels <- rep(c("(Intercept)", attr(terms(fMod), "term.labels")[attr(terms(fMod), "term.labels") %in% tt$between]), tableTerms)
    term <- c("(Intercept)", attr(terms(fMod), "term.labels")[attr(terms(fMod), "term.labels") %in% tt$between])
    rowIndex <- list()
    for(i in 1:length(term)){
      rowIndex[[i]] <- which(termsLabels == term[i])
      names(rowIndex)[i] <- term[i]
    }
    Lbetween <- lapply(rowIndex, function(x){
      Lbetween[x,, drop = FALSE]
    })
  }
  Qbetween <- Qwithin <- sapply(Lbetween, function(x){
    Q <- Matrix::rankMatrix(x %*% MASS::ginv(t(X) %*% X) %*% t(x))
    attributes(Q) <- NULL
    Q
  })

  if(waldF == FALSE){
    sBetween <- apply(cbind(pBetween, Qbetween), 1, min)
    X2Between <- v * UBetween
    X2Within <- v * UWithin
    Qwithin <- Qwithin[1:length(tt$within)]
  }else{
    if(length(tt$between) == 0L){
      Qbetween <- pBetween <- sBetween <- NULL
    }else{
      pqBetween <- waldBetween$NumDF
      Qbetween <- Qbetween[-1]
      pBetween <- pqBetween / Qbetween
      #names(pBetween) <- tt$between
      sBetween <- apply(cbind(pBetween, Qbetween), 1, min)
    }
    pqWithin <- waldWithin$NumDF
    Qwithin <- Qwithin[1:length(tt$within)]
    pWithin <- pqWithin / Qwithin
    names(pWithin) <- tt$within
  }
  sWithin <- apply(cbind(pWithin, Qwithin), 1, min)[1:length(tt$within)]

  U <- list(EB = E$EB, HB = H$HB, UBetween = UBetween,
            sBetween = sBetween, pBetween = pBetween, qBetween = Qbetween,
            X2Between = X2Between, E = E$E, H = H$H,
            UWithin = UWithin, sWithin = sWithin,
            pWithin = pWithin, qWithin = Qwithin, X2Within = X2Within,
            X = X, v = v, Lbetween = Lbetween, Lwithin = Lwithin, L = L, n = n, type = type, waldF = waldF)

  class(U) <- c("list", "U.ulm")
  return(U)
}

#' Multivariate Tests: Hotelling-Lawley Trace
#'
#' @export
#'
#' @description
#' This function will calculate Hotelling-Lawley Trace \eqn{U} statistics from marginal models,
#' mixed-effects models, and multivariate linear models.
#'
#' @details
#'
#' **H and E matrices**
#'
#' The Hotelling-Lawley trace statistic for RM ANOVA is calculated as
#'
#' \eqn{U = trace(}**E**\eqn{ ^ {-1}}**H**\eqn{)}
#'
#' where the SSCP error matrix is
#'
#' **E = M'** \eqn{(}**Y'Y**\eqn{-} **\eqn{\hat{B}}'** \eqn{(}**X'X**\eqn{)}**\eqn{\hat{B}}** \eqn{)}**M**
#'
#' and the SSCP hypothesis matrix is
#'
#' **H** \eqn{=} **M'** \eqn{(}**L\eqn{\hat{B}}** \eqn{)}'\eqn{(}**L**\eqn{(}**X'X**\eqn{) ^ {-}}**L'** \eqn{) ^ {-1}}**L\eqn{\hat{B}}M**
#'
#' where **M** is an identity matrix for between-subject effects and a
#' sum-to-zero contrast matrix for within-subject effects, **Y** is a combined
#' matrix of observations, **\eqn{\hat{B}}** is a matrix of combined coefficients from the
#' ordinary least squares models, and **X** is the design matrix of the
#' independent variables.
#'
#' For marginal and mixed effects models, the **E** matrix is calculated from
#' the marginal variance-covariance matrix of the random/repeated effects
#' (**V**)
#'
#' **E =** \eqn{n}**M'VM**
#'
#' where **V = ZGZ' + R** (**Z** is the random effects design matrix, **G** is
#' the corresponding (co)variances of the random effects matrix, and **R**
#' represents the (co)variances of the repeated effects matrix).  The **H**
#' matrix is calculated with the same formula as above.  The **\eqn{\hat{B}}** matrix is
#' calculated by switching the reference levels of the within-subject effect,
#' refitting the models, and extracting the model coefficients.
#'
#' **Wald**
#'
#' For marginal and mixed-effects models, we start with a Wald-type quadratic
#' form
#'
#' \eqn{F = (}**\eqn{\hat{\beta}}'L**\eqn{(}**L'** \eqn{(}**X'V**\eqn{ ^ {-1}}**X**\eqn{) ^ {-}}**L**\eqn{) ^ {-1}}**L'\eqn{\hat{\beta}}** \eqn{) / rank(}**L'** \eqn{(}**X'V**\eqn{ ^ {-1}}**X**\eqn{) ^ {-}}**L**\eqn{)}
#'
#' where **\eqn{\hat{\beta}}** is a vector of fixed effects coefficients, **X** is the fixed effects design matrix, and
#' **V** is the marginal variance-covariance matrix
#' where **V = ZGZ' + R** (**Z** is the random effects design matrix, **G** is
#' the corresponding (co)variances of the random effects matrix, and **R**
#' represents the (co)variances of the repeated effects matrix).  The contrast coefficients matrix (**L**) is calculated from the **LU**
#' decomposition of the crossproducts matrix (**X'X**).
#'
#' It follows that the Wald-\eqn{\chi ^ 2} statistic is the numerator degrees of
#' freedom (\eqn{rank(}**L'** \eqn{(}**X'V**\eqn{ ^ {-1}}**X**\eqn{) ^ {-}}**L**\eqn{)})
#' multiplied by the F-statistic.  The Hotelling-Lawley trace (\eqn{U})
#' is then calculated as
#'
#' \eqn{U = \chi ^ 2 / (n - rank(}**X_b** \eqn{))}.
#'
#' Let \eqn{v} be the residual degrees of freedom, \eqn{p = rank(}**H**\eqn{+}**E**\eqn{)},
#' \eqn{q = rank(}**L**\eqn{(}**X'X**\eqn{)^{-1}}**L'** \eqn{)}, \eqn{s = min(p, q)}, \eqn{m = (|p - q| - 1) / 2} and \eqn{n = (v - p - 1) / 2}.
#'
#' There are two approximations used for marginal and mixed-effects models:
#' Pillai-Samson (1959) and McKeon (1974).  The Pillai-Samson approximation is recommended for
#' \eqn{n \le 0} while the McKeon approximation is only defined for \eqn{n > 0}.
#' The McKeon approximation is said to be more accurate overall.
#'
#' **McKeon (1974)**
#'
#' The McKeon statistic approximately follows an F-distribution where
#' \deqn{F_U = (U / c)((4 + (pq + 2) / (b - 1)) / (pq))}
#' with \eqn{pq} numerator degrees of freedom and \eqn{4 + (pq + 2) / (b - 1)}
#' denominator degrees of freedom, where
#'
#' \eqn{b = (p + 2n)(q + 2n) / (2(2n + 1)(n - 1))} and
#'
#' \eqn{c = (2 + (pq + 2) / (b - 1)) / (2n)}
#'
#' If \eqn{s = 1}, the F-statistic is exact.
#'
#' **Pillai-Samson (1959)**
#'
#' The Pillai-Samson statistic approximately follows an F-distribution where
#' \deqn{F_U = (2(sn + 1)U) / (s ^ 2(2m + s + 1))} with \eqn{s(2m + s + 1)}
#' numerator degrees of freedom and \eqn{2(sn + 1)} denominator degrees of
#' freedom.
#'
#' If \eqn{s = 1}, the F-statistic is exact.
#'
#' @param fMod A model of class \code{\link[nlme]{gls}}, \code{\link[nlme]{lme}}, or \code{mlm}.
#' @param individual The subject whose marginal variance-covariance matrix
#' should be used, where **V = ZGZ' + R** (**Z** is the random effects design
#' matrix, **G** is the corresponding (co)variances of the random effects matrix,
#' and **R** represents the (co)variances of the repeated effects matrix).  The
#' default is the first subject in the factor.
#' @param approximation Approximation used for Hotelling-Lawley trace.  Options are
#' McKeon F ("McKeon"), Pillai-Samson F ("Pillai-Samson"), and Wald \eqn{\chi ^ 2}
#' ("Wald").
#' @param within The name that should be used for the within-subject effect.
#' The default is "Time".
#' @param waldF Wald F-tests are used to calculate the Wald \eqn{\chi ^ 2} and
#' Hotelling-Lawley trace statistics instead of the **H** and **E** matrices.
#'
#' @references {\url{https://support.sas.com/resources/papers/proceedings/proceedings/sugi23/Stats/p229.pdf}}
#' @references {\url{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.4/statug/statug_introreg_sect038.htm}}
#' @references {McKeon, J. J. (1974). "F Approximations to the Distribution of Hotelling's \eqn{T_0^2}." Biometrika 61:381-383.}
#' @references {Pillai, K. C. S., and Samson, P., Jr. (1959). "On Hotelling's Generalization of \eqn{T_0^2}." Biometrika 46:160-168.}
#' @references {Wright, S. P. (1994). Adjusted F Tests for Repeated Measures with
#' the MIXED Procedure. Knoxville: Statistics Department, University of Tennessee.
#' \url{https://support.sas.com/resources/papers/proceedings-archive/SUGI95/Sugi-95-195%20Wright.pdf}}
#' @references {\url{https://brainder.org/tag/lawley-hotelling-trace/}}
#'
#' @returns
#' This function will return two lists: a within-subject effects list
#' (`withinTests`) and a between-subject effects list (`betweenTests`) if
#' between-subject effects are in the model.
#' @returns `withinTests` is a list containing: \tabular{ll}{
#'    `hlt`  \tab A `data.frame` containing the Hotelling-Lawley Trace \eqn{U} statistics and
#'    hypothesis tests for these statistics. \cr
#'    `parmsHL`  \tab A `data.frame` containing the values used to calculate
#'    the Hotelling-Lawley Trace \eqn{U} statistics. \cr
#'    `approximation`  \tab Approximation used for the Hotelling-Lawley Trace \eqn{U} statistics.  Options are McKeon F ("McKeon"),
#' Pillai-Samson F ("Pillai-Samson"), and Wald \eqn{\chi ^ 2} ("Wald"). \cr
#'    `E`  \tab The sums of squares and crossproducts error matrix for within-subject effects. \cr
#'    `H`  \tab The sums of squares and crossproducts hypothesis matrices for within-subject effects. \cr
#'    `L`  \tab The contrast coefficients matrices. \cr
#'    `type`  \tab The type of **L** matrix used.  Type 1 is sequential. \cr
#' \tab \cr
#' }
#' `betweenTests` is a list containing: \tabular{ll}{
#'    `hlt`  \tab A `data.frame` containing the Hotelling-Lawley Trace \eqn{U} statistics and
#'    hypothesis tests for these statistics. \cr
#'    `parmsHL`  \tab A `data.frame` containing the values used to calculate
#'    the Hotelling-Lawley Trace \eqn{U} statistics. \cr
#'    `approximation`  \tab Approximation used for the Hotelling-Lawley Trace \eqn{U} statistics.  Options are McKeon F ("McKeon"),
#' Pillai-Samson F ("Pillai-Samson"), and Wald \eqn{\chi ^ 2} ("Wald"). \cr
#'    `EB`  \tab The sums of squares and crossproducts error matrix for between-subject effects. \cr
#'    `HB`  \tab The sums of squares and crossproducts hypothesis matrices for between-subject effects. \cr
#'    `L`  \tab The contrast coefficients matrices. \cr
#'    `type`  \tab The type of **L** matrix used.  Type 1 is sequential. \cr
#'    `Vmatrix`  \tab The marginal variance-covariance matrix of the random/repeated effects. \cr
#' \tab \cr
#' }
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#' # Chi-square statistic
#' hlTrace(glsMod, approximation = "Wald")
#' # Pillai-Samson F statistic
#' hlTrace(glsMod, approximation = "P")
#' # McKeon F statistic
#' hlTrace(glsMod)
#'
#' # Linear mixed-effects model
#' # Unstructured covariance
#' lmeMod <- lme(Distance ~ Sex * Age,
#'               random = ~ Age | Subject,
#'               na.action = na.omit,
#'               data = orthodontLong)
#' # Chi-square statistic
#' hlTrace(lmeMod, approximation = "Wald")
#' # Pillai-Samson F statistic
#' hlTrace(lmeMod, approximation = "P")
#' # McKeon F statistic
#' hlTrace(lmeMod)
#'
#' # Multivariate linear model
#' manovaMod <- manova(cbind(Distance8, Distance10, Distance12, Distance14) ~ Sex,
#'                     data = orthodontWide)
#' # McKeon F statistic
#' hlTrace(manovaMod)
#' # Pillai-Samson F statistic
#' hlTrace(manovaMod, approximation = "P")
#' # Chi-square statistic
#' hlTrace(manovaMod, approximation = "Wald")
#'
#' @seealso
#' {\code{\link[mmemmuris]{Hmatrix}}, \code{\link[mmemmuris]{coefs}}, \code{\link[mmemmuris]{Vmatrix}}, \code{\link[mmemmuris]{Ematrix}}, \code{\link[mmemmuris]{waldF}}, \code{\link[mmemmuris]{MLtoREML}}}

hlTrace <- function(fMod, individual = NULL, approximation = c("McKeon", "Pillai-Samson", "Wald"),
                        within = "Time", waldF = TRUE, ss = NULL){
  approximation <- match.arg(approximation)
  if(any(class(fMod) %in% c("gls", "lme"))){
    hlt <- mmemmuris:::hlTrace.ulm(fMod, individual, approximation, waldF = waldF, ss = ss)
    return(hlt)
  }else if(any(class(fMod) %in% "mlm")){
    hlt <- mmemmuris:::hlTrace.mlm(fMod, approximation, within)
    return(hlt)
  }else
    stop('Please provide a model of class "gls", "lme", or "mlm".', call. = FALSE)
}

hlTrace.ulm <- function(fMod, individual = NULL, approximation = c("McKeon", "Pillai-Samson", "Wald"), waldF = TRUE, ss = NULL){
  approximation <- match.arg(approximation)
  tt <- mmemmuris::termsType(fMod)
  hlt <- mmemmuris:::U.ulm(fMod, individual, waldF = waldF, ss = ss)

  if(approximation == "McKeon"){
    mckeonFApprox <- mmemmuris:::mckeonF(hlt)

    withinTests <- mckeonFApprox$hltWithin

    if(!is.null(mckeonFApprox$hltBetween)){
      betweenTests <- mckeonFApprox$hltBetween

      hlt <- list(betweenTests = betweenTests,
                  withinTests = withinTests,
                  waldF = hlt$waldF)
      class(hlt) <- c("list", "multiTest", "hlt", "ulm")
    }else{
      hlt <- list(withinTests = withinTests,
                  waldF = hlt$waldF)
      class(hlt) <- c("list", "multiTest", "hlt", "ulm")
    }
    return(hlt)
  }else if(approximation == "Pillai-Samson"){
    psFApprox <- mmemmuris:::psF(hlt)

    withinTests <- psFApprox$hltWithin

    if(!is.null(psFApprox$hltBetween)){
      betweenTests <- psFApprox$hltBetween

      hlt <- list(betweenTests = betweenTests,
                  withinTests = withinTests,
                  waldF = hlt$waldF)
      class(hlt) <- c("list", "multiTest", "hlt", "ulm")
    }else{
      hlt <- list(withinTests = withinTests,
                  waldF = hlt$waldF)
      class(hlt) <- c("list", "multiTest", "hlt", "ulm")
    }
    return(hlt)
  }else if(approximation == "Wald"){
    waldChiApprox <- mmemmuris:::waldChi(hlt)

    withinTests <- waldChiApprox$hltWithin

    if(!is.null(waldChiApprox$hltBetween)){
      betweenTests <- waldChiApprox$hltBetween

      hlt <- list(betweenTests = betweenTests,
                  withinTests = withinTests,
                  waldF = hlt$waldF)
      class(hlt) <- c("list", "multiTest", "hlt", "ulm")
    }else{
      hlt <- list(withinTests = withinTests,
                  waldF = hlt$waldF)
      class(hlt) <- c("list", "multiTest", "hlt", "ulm")
    }
    return(hlt)
  }
}

hlTrace.mlm <- function(fMod, approximation = c("McKeon", "Pillai-Samson", "Wald"),
                        within = "Time"){
  approximation <- match.arg(approximation)
  hlt <- mmemmuris:::U.mlm(fMod, within)

  if(approximation == "McKeon"){
    mckeonFApprox <- mmemmuris:::mckeonF(hlt)

    withinTests <- mckeonFApprox$hltWithin

    if(!is.null(mckeonFApprox$hltBetween)){
      betweenTests <- mckeonFApprox$hltBetween
      hlt <- list(betweenTests = betweenTests,
                  withinTests = withinTests,
                  waldF = FALSE)
      class(hlt) <- c("list", "multiTest", "hlt", "mlm")
    }else{
      hlt <- list(withinTests = withinTests,
                  waldF = FALSE)
      class(hlt) <- c("list", "multiTest", "hlt", "mlm")
    }
    return(hlt)
  }else if(approximation == "Pillai-Samson"){
    psFApprox <- mmemmuris:::psF(hlt)

    withinTests <- psFApprox$hltWithin

    if(!is.null(psFApprox$hltBetween)){
      betweenTests <- psFApprox$hltBetween
      hlt <- list(betweenTests = betweenTests,
                  withinTests = withinTests,
                  waldF = FALSE)
      class(hlt) <- c("list", "multiTest", "hlt", "mlm")
    }else{
      hlt <- list(withinTests = withinTests,
                  waldF = FALSE)
      class(hlt) <- c("list", "multiTest", "hlt", "mlm")
    }
    return(hlt)
  }else if(approximation == "Wald"){
    waldChiApprox <- mmemmuris:::waldChi(hlt)

    withinTests <- waldChiApprox$hltWithin

    if(!is.null(waldChiApprox$hltBetween)){
      betweenTests <- waldChiApprox$hltBetween
      hlt <- list(betweenTests = betweenTests,
                  withinTests = withinTests,
                  waldF = FALSE)
      class(hlt) <- c("list", "multiTest", "hlt", "mlm")
    }else{
      hlt <- list(withinTests = withinTests,
                  waldF = FALSE)
      class(hlt) <- c("list", "multiTest", "hlt", "mlm")
    }
    return(hlt)
  }
}

makeHLTrace.ulm <- function(U, approximation = c("McKeon", "Pillai-Samson", "Wald")){
  approximation <- match.arg(approximation)

  if(approximation == "McKeon"){
    mckeonFApprox <- mckeonF(U)

    withinTests <- mckeonFApprox$hltWithin

    if(!is.null(mckeonFApprox$hltBetween)){
      betweenTests <- mckeonFApprox$hltBetween

      hlt <- list(betweenTests = betweenTests,
                  withinTests = withinTests,
                  waldF = U$waldF)
      class(hlt) <- c("list", "multiTest", "hlt", "ulm")
    }else{
      hlt <- list(withinTests = withinTests,
                  waldF = U$waldF)
      class(hlt) <- c("list", "multiTest", "hlt", "ulm")
    }
    return(hlt)
  }else if(approximation == "Pillai-Samson"){
    psFApprox <- psF(U)

    withinTests <- psFApprox$hltWithin

    if(!is.null(psFApprox$hltBetween)){
      betweenTests <- psFApprox$hltBetween

      hlt <- list(betweenTests = betweenTests,
                  withinTests = withinTests,
                  waldF = U$waldF)
      class(hlt) <- c("list", "multiTest", "hlt", "ulm")
    }else{
      hlt <- list(withinTests = withinTests,
                  waldF = U$waldF)
      class(hlt) <- c("list", "multiTest", "hlt", "ulm")
    }
    return(hlt)
  }else if(approximation == "Wald"){
    waldChiApprox <- waldChi(U)

    withinTests <- waldChiApprox$hltWithin

    if(!is.null(waldChiApprox$hltBetween)){
      betweenTests <- waldChiApprox$hltBetween

      hlt <- list(betweenTests = betweenTests,
                  withinTests = withinTests,
                  waldF = U$waldF)
      class(hlt) <- c("list", "multiTest", "hlt", "ulm")
    }else{
      hlt <- list(withinTests = withinTests,
                  waldF = U$waldF)
      class(hlt) <- c("list", "multiTest", "hlt", "ulm")
    }
    return(hlt)
  }
}

#' Multivariate Tests: Roy's Greatest Root
#'
#' @export
#'
#' @description
#' This function will calculate Roy's Greatest Root \eqn{\Theta} statistics from marginal models,
#' mixed-effects models, and multivariate linear models.
#'
#' @details
#'
#' The Roy's greatest root statistic for RM ANOVA is calculated as
#'
#' \eqn{\Theta = max(\lambda_i)} of **E**\eqn{ ^ {-1}}**H**
#'
#' where the SSCP error matrix is
#'
#' **E = M'** \eqn{(}**Y'Y**\eqn{-} **\eqn{\hat{B}}'** \eqn{(}**X'X**\eqn{)}**\eqn{\hat{B}}** \eqn{)}**M**
#'
#' and the SSCP hypothesis matrix is
#'
#' **H** \eqn{=} **M'** \eqn{(}**L\eqn{\hat{B}}** \eqn{)}'\eqn{(}**L**\eqn{(}**X'X**\eqn{) ^ {-}}**L'** \eqn{) ^ {-1}}**L\eqn{\hat{B}}M**
#'
#' where **M** is an identity matrix for between-subject effects and a
#' sum-to-zero contrast matrix for within-subject effects, **Y** is a combined
#' matrix of observations, **\eqn{\hat{B}}** is a matrix of combined coefficients from the
#' ordinary least squares models, and **X** is the design matrix of the
#' independent variables.
#'
#' For marginal and mixed effects models, the **E** matrix is calculated from
#' the marginal variance-covariance matrix of the random/repeated effects
#' (**V**)
#'
#' **E =** \eqn{n}**M'VM**
#'
#' where **V = ZGZ' + R** (**Z** is the random effects design matrix, **G** is
#' the corresponding (co)variances of the random effects matrix, and **R**
#' represents the (co)variances of the repeated effects matrix).  The **H**
#' matrix is calculated with the same formula as above.  The **\eqn{\hat{B}}** matrix is
#' calculated by switching the reference levels of the within-subject effect,
#' refitting the models, and extracting the model coefficients.
#'
#' **Wald**
#'
#' For marginal and mixed-effects models, we start with a Wald-type quadratic
#' form
#'
#' \eqn{F = (}**\eqn{\hat{\beta}}'L**\eqn{(}**L'** \eqn{(}**X'V**\eqn{ ^ {-1}}**X**\eqn{) ^ {-}}**L**\eqn{) ^ {-1}}**L'\eqn{\hat{\beta}}** \eqn{) / rank(}**L'** \eqn{(}**X'V**\eqn{ ^ {-1}}**X**\eqn{) ^ {-}}**L**\eqn{)}
#'
#' where **\eqn{\hat{\beta}}** is a vector of fixed effects coefficients, **X** is the fixed effects design matrix, and
#' **V** is the marginal variance-covariance matrix
#' where **V = ZGZ' + R** (**Z** is the random effects design matrix, **G** is
#' the corresponding (co)variances of the random effects matrix, and **R**
#' represents the (co)variances of the repeated effects matrix).  The contrast coefficients matrix (**L**) is calculated from the **LU**
#' decomposition of the crossproducts matrix (**X'X**).
#'
#' It follows that the Wald-\eqn{\chi ^ 2} statistic is the numerator degrees of
#' freedom (\eqn{rank(}**L'** \eqn{(}**X'V**\eqn{ ^ {-1}}**X**\eqn{) ^ {-}}**L**\eqn{)})
#' multiplied by the F-statistic.  The Hotelling-Lawley trace (\eqn{U})
#' is then calculated as
#'
#' \eqn{U = \chi ^ 2 / (n - rank(}**X_b** \eqn{))}.
#'
#' For marginal and mixed-effects models, if \eqn{s = min(p, q) = 1}, Roy's
#' greatest root is calculated as \eqn{\Theta = U}.  With `waldF = TRUE`,
#' Roy's greatest root is only calculated when \eqn{s = 1} for within-subject
#' effects.
#'
#' Let \eqn{v} be the residual degrees of freedom, \eqn{p = rank(}**H**\eqn{+}**E**\eqn{)},
#' \eqn{q = rank(}**L**\eqn{(}**X'X**\eqn{)^{-1}}**L'** \eqn{)}, and \eqn{s = min(p, q)}.
#'
#' The statistic approximately follows an F-distribution
#' \deqn{F_\Theta = \Theta((v - max(p, q) + q) / max(p, q))}
#' with \eqn{max(p, q)} numerator degrees of freedom and \eqn{v - max(p, q) + q}
#' denominator degrees of freedom.
#'
#' @param fMod A model of class \code{\link[nlme]{gls}} or \code{\link[nlme]{lme}}.
#' @param individual The subject whose marginal variance-covariance matrix
#' should be used, where **V = ZGZ' + R** (**Z** is the random effects design
#' matrix, **G** is the corresponding (co)variances of the random effects matrix,
#' and **R** represents the (co)variances of the repeated effects matrix).  The
#' default is the first subject in the factor.
#' @param within The name that should be used for the within-subject effect.
#' The default is "Time".
#' @param waldF Wald F-tests are used to calculate the Wald \eqn{\chi ^ 2} and
#' Hotelling-Lawley trace statistics instead of the **H** and **E** matrices.
#'
#' @returns
#' This function will return two lists: a within-subject effects list
#' (`withinTests`) and a between-subject effects list (`betweenTests`) if
#' between-subject effects are in the model.
#' @returns `withinTests` is a list containing: \tabular{ll}{
#'    `rgr`  \tab A `data.frame` containing the Roy's Greatest Root \eqn{\Theta} statistics and
#'    hypothesis tests for these statistics. \cr
#'    `parmsRGR`  \tab A `data.frame` containing the values used to calculate
#'    the Roy's Greatest Root \eqn{\Theta} statistics. \cr
#'    `E`  \tab The sums of squares and crossproducts error matrix for within-subject effects. \cr
#'    `H`  \tab The sums of squares and crossproducts hypothesis matrices for within-subject effects. \cr
#'    `L`  \tab The contrast coefficients matrices. \cr
#'    `type`  \tab The type of **L** matrix used.  Type 1 is sequential. \cr
#' \tab \cr
#' }
#' `betweenTests` is a list containing: \tabular{ll}{
#'    `rgr`  \tab A `data.frame` containing the Roy's Greatest Root \eqn{\Theta} statistics and
#'    hypothesis tests for these statistics. \cr
#'    `parmsRGR`  \tab A `data.frame` containing the values used to calculate
#'    the Roy's Greatest Root \eqn{\Theta} statistics. \cr
#'    `EB`  \tab The sums of squares and crossproducts error matrix for between-subject effects. \cr
#'    `HB`  \tab The sums of squares and crossproducts hypothesis matrices for between-subject effects. \cr
#'    `L`  \tab The contrast coefficients matrices. \cr
#'    `type`  \tab The type of **L** matrix used.  Type 1 is sequential. \cr
#'    `Vmatrix`  \tab The marginal variance-covariance matrix of the random/repeated effects. \cr
#' \tab \cr
#' }
#'
#' @references {\url{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.4/statug/statug_introreg_sect038.htm}}
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#' # F statistic
#' royGR(glsMod)
#'
#' # Linear mixed-effects model
#' # Unstructured covariance
#' lmeMod <- lme(Distance ~ Sex * Age,
#'               random = ~ Age | Subject,
#'               na.action = na.omit,
#'               data = orthodontLong)
#' # F statistic
#' royGR(lmeMod)
#'
#' # Multivariate linear model
#' manovaMod <- manova(cbind(Distance8, Distance10, Distance12, Distance14) ~ Sex,
#'                     data = orthodontWide)
#' # F statistic
#' royGR(manovaMod)
#'
#' @seealso
#' {\code{\link[mmemmuris]{Hmatrix}}, \code{\link[mmemmuris]{coefs}}, \code{\link[mmemmuris]{Vmatrix}}, \code{\link[mmemmuris]{Ematrix}}, \code{\link[mmemmuris]{waldF}}, \code{\link[mmemmuris]{MLtoREML}}}

royGR <- function(fMod, individual = NULL, within = "Time", waldF = TRUE, ss = NULL){
  if(any(class(fMod) %in% c("gls", "lme"))){
    rgr <- mmemmuris:::royGR.ulm(fMod, individual, waldF = waldF, ss = ss)
    return(rgr)
  }else if(any(class(fMod) %in% "mlm")){
    rgr <- mmemmuris:::royGR.mlm(fMod, within)
    return(rgr)
  }else
    stop('Please provide a model of class "gls", "lme", or "mlm".', call. = FALSE)
}

royGR.mlm <- function(fMod, within = "Time"){
  theta <- mmemmuris:::theta.mlm(fMod, within)

  rWithin <- apply(cbind(theta$pWithin, theta$qWithin), 1, max)
  parmsRGRWithin <- data.frame(p = theta$pWithin, q = theta$qWithin, s = theta$sWithin,
                               r = rWithin, v = theta$v)
  FrgrWithin <- theta$thetaWithin * ((theta$v - rWithin + theta$qWithin) / rWithin)
  termsRGRWithin <- data.frame(dfNum = rWithin,
                               dfDen = theta$v - rWithin + theta$qWithin,
                               theta = theta$thetaWithin, F = FrgrWithin,
                               p.value = 1 - pf(FrgrWithin, rWithin, theta$v - rWithin + theta$qWithin))

  colnames(termsRGRWithin) <- c("Num df", "Den df", "theta", "F", "Pr(>F)")
  termsRGRWithin <- cbind(termsRGRWithin, mmemmuris::sigStars(termsRGRWithin$`Pr(>F)`))
  colnames(termsRGRWithin)[length(colnames(termsRGRWithin))] <- ""

  withinTests <- list(rgr = termsRGRWithin,
                      parmsRGR = parmsRGRWithin,
                      approximation = "Pillai",
                      E = theta$E, H = theta$H,
                      Lwithin = theta$Lwithin, type = theta$type)
  if(length(theta$thetaBetween[names(theta$thetaBetween) != "(Intercept)"]) > 0L){
    rBetween <- apply(cbind(theta$pBetween, theta$qBetween), 1, max)
    parmsRGRBetween <- data.frame(p = theta$pBetween, q = theta$qBetween,
                                  s = theta$sBetween, r = rBetween,
                                  v = theta$v)
    FrgrBetween <- theta$thetaBetween * ((theta$v - rBetween + theta$qBetween) / rBetween)
    termsRGRBetween <- data.frame(dfNum = rBetween,
                                  dfDen = theta$v - rBetween + theta$qBetween,
                                  theta = theta$thetaBetween, F = FrgrBetween,
                                  p.value = 1 - pf(FrgrBetween, rBetween, theta$v - rBetween + theta$qBetween))

    colnames(termsRGRBetween) <- c("Num df", "Den df", "theta", "F", "Pr(>F)")
    termsRGRBetween <- cbind(termsRGRBetween, mmemmuris::sigStars(termsRGRBetween$`Pr(>F)`))
    colnames(termsRGRBetween)[length(colnames(termsRGRBetween))] <- ""

    betweenTests <- list(rgr = termsRGRBetween[rownames(termsRGRBetween) != "(Intercept)", ],
                         parmsRGR = parmsRGRBetween[rownames(parmsRGRBetween) != "(Intercept)", ],
                         approximation = "Pillai",
                         EB = theta$EB, HB = theta$HB,
                         Lbetween = theta$Lbetween, L = theta$L, type = theta$type,
                         Vmatrix = theta$EB / theta$n)

    rgr <- list(betweenTests = betweenTests,
                withinTests = withinTests,
                waldF = FALSE)
    class(rgr) <- c("list", "multiTest", "rgr", "mlm")
  }else{
    rgr <- list(withinTests = withinTests,
                waldF = FALSE)
    class(rgr) <- c("list", "multiTest", "rgr", "mlm")
  }
  return(rgr)
}

royGR.ulm <- function(fMod, individual = NULL, waldF = TRUE, ss = NULL){
  tt <- mmemmuris::termsType(fMod)
  theta <- mmemmuris:::theta.ulm(fMod, individual, waldF = waldF, ss = ss)

  rWithin <- apply(cbind(theta$pWithin, theta$qWithin), 1, max)
  parmsRGRWithin <- data.frame(p = theta$pWithin, q = theta$qWithin, s = theta$sWithin,
                               r = rWithin, v = theta$v)
  FrgrWithin <- theta$thetaWithin * ((theta$v - rWithin + theta$qWithin) / rWithin)
  termsRGRWithin <- data.frame(dfNum = rWithin,
                               dfDen = theta$v - rWithin + theta$qWithin,
                               theta = theta$thetaWithin, F = FrgrWithin,
                               p.value = 1 - pf(FrgrWithin, rWithin, theta$v - rWithin + theta$qWithin))

  colnames(termsRGRWithin) <- c("Num df", "Den df", "theta", "F", "Pr(>F)")
  termsRGRWithin <- cbind(termsRGRWithin, mmemmuris::sigStars(termsRGRWithin$`Pr(>F)`))
  colnames(termsRGRWithin)[length(colnames(termsRGRWithin))] <- ""

  if(theta$waldF == TRUE){
    withinTests <- list(rgr = termsRGRWithin[theta$sWithin == 1, ],
                        parmsRGR = parmsRGRWithin,
                        approximation = "Pillai",
                        E = theta$E, H = theta$H,
                        Lwithin = theta$Lwithin, type = theta$type)
  }else{
    withinTests <- list(rgr = termsRGRWithin,
                        parmsRGR = parmsRGRWithin,
                        approximation = "Pillai",
                        E = theta$E, H = theta$H,
                        Lwithin = theta$Lwithin, type = theta$type)
  }

  if(!is.null(theta$thetaBetween)){
    rBetween <- apply(cbind(theta$pBetween, theta$qBetween), 1, max)
    parmsRGRBetween <- data.frame(p = theta$pBetween, q = theta$qBetween,
                                  s = theta$sBetween, r = rBetween,
                                  v = theta$v)
    FrgrBetween <- theta$thetaBetween * ((theta$v - rBetween + theta$qBetween) / rBetween)
    termsRGRBetween <- data.frame(dfNum = rBetween,
                                  dfDen = theta$v - rBetween + theta$qBetween,
                                  theta = theta$thetaBetween, F = FrgrBetween,
                                  p.value = 1 - pf(FrgrBetween, rBetween, theta$v - rBetween + theta$qBetween))

    colnames(termsRGRBetween) <- c("Num df", "Den df", "theta", "F", "Pr(>F)")
    termsRGRBetween <- cbind(termsRGRBetween, mmemmuris::sigStars(termsRGRBetween$`Pr(>F)`))
    colnames(termsRGRBetween)[length(colnames(termsRGRBetween))] <- ""

    if(theta$waldF){
      betweenTests <- list(rgr = termsRGRBetween[rownames(termsRGRBetween) != "(Intercept)" & theta$sBetween == 1, ],
                           parmsRGR = parmsRGRBetween[rownames(parmsRGRBetween) != "(Intercept)", ],
                           approximation = "Pillai",
                           EB = theta$EB, HB = theta$HB,
                           Lbetween = theta$Lbetween, L = theta$L, type = theta$type,
                           Vmatrix = theta$EB / theta$n)
    }else{
      betweenTests <- list(rgr = termsRGRBetween[rownames(termsRGRBetween) != "(Intercept)", ],
                           parmsRGR = parmsRGRBetween[rownames(parmsRGRBetween) != "(Intercept)", ],
                           approximation = "Pillai",
                           EB = theta$EB, HB = theta$HB,
                           Lbetween = theta$Lbetween, L = theta$L, type = theta$type,
                           Vmatrix = theta$EB / theta$n)
    }

    if(nrow(betweenTests$rgr) > 0L)
      rgr <- list(betweenTests = betweenTests, withinTests = withinTests,
                  waldF = theta$waldF)
    else
      rgr <- list(withinTests = withinTests,
                  waldF = theta$waldF)
    class(rgr) <- c("list", "multiTest", "rgr", "ulm")
  }else{
    rgr <- list(withinTests = withinTests,
                waldF = theta$waldF)
    class(rgr) <- c("list", "multiTest", "rgr", "ulm")
  }
  return(rgr)
}

theta.mlm <- function(fMod, within = "Time"){
  type <- "1"
  E <- mmemmuris:::Ematrix.mlm(fMod)
  n <- E$n
  H <- mmemmuris:::Hmatrix.mlm(fMod, within)

  # Calculate Roy's greatest root statistic
  thetaBetween <- sapply(H$HB, function(x){
    max(Re(eigen(solve(E$EB) %*% x)$values))
  })
  thetaWithin <- sapply(H$H, function(x){
    max(Re(eigen(solve(E$E) %*% x)$values))
  })

  pBetween <- sapply(H$HB, function(x){
    pBetween <- Matrix::rankMatrix(x + E$EB)
    attributes(pBetween) <- NULL
    pBetween
  })
  pWithin <- sapply(H$H, function(x){
    pWithin <- Matrix::rankMatrix(x + E$E)
    attributes(pWithin) <- NULL
    pWithin
  })

  X <- H$X
  Qbetween <- Qwithin <- sapply(H$L, function(x){
    Q <- Matrix::rankMatrix(x %*% MASS::ginv(t(X) %*% X) %*% t(x))
    attributes(Q) <- NULL
    Q
  })

  rankX <- Matrix::rankMatrix(X)
  attributes(rankX) <- NULL
  v <- n - rankX

  sBetween <- apply(cbind(pBetween, Qbetween), 1, min)
  sWithin <- apply(cbind(pWithin, Qwithin), 1, min)

  theta <- list(EB = E$EB, HB = H$HB, thetaBetween = thetaBetween,
                sBetween = sBetween, pBetween = pBetween, qBetween = Qbetween,
                E = E$E, H = H$H,
                thetaWithin = thetaWithin, sWithin = sWithin,
                pWithin = pWithin, qWithin = Qwithin, X = X, v = v,
                Lbetween = H$L, n = n, type = type)

  class(theta) <- c("list", "theta.mlm")
  return(theta)
}

theta.ulm <- function(fMod, individual = NULL, waldF = TRUE, ss = NULL){
  type <- "1"
  if(mmemmuris::covStruct(fMod) %in% c("cs", "other"))
    stop("Only unstructured covariance models are allowed.", call. = FALSE)
  tt <- mmemmuris::termsType(fMod)
  EBsingular <- Esingular <- FALSE
  E <- mmemmuris:::Ematrix.ulm(fMod, individual, ss)
  H <- mmemmuris:::Hmatrix.ulm(fMod)
  n <- E$n
  if("error" %in% class(tryCatch(any(sapply(H$HB, function(x){ Matrix::rankMatrix(x + E$EB) })), error = function(e) { e })))
    stop('Please pick a different E matrix with the "individual" argument.', call. = FALSE)
  if("error" %in% class(tryCatch(any(sapply(H$H, function(x){ Matrix::rankMatrix(x + E$E) })), error = function(e) { e })))
    stop('Please pick a different E matrix with the "individual" argument.', call. = FALSE)
  if(is.null(mmemmuris:::inverseMatrix(E$EB)))
    EBsingular <- TRUE
  if(is.null(mmemmuris:::inverseMatrix(E$E)))
    Esingular <- TRUE
  if(EBsingular == TRUE | Esingular == TRUE){
    waldF <- TRUE
    warning("EB ^ {-1} or E ^ {-1} are singular matrices.  Wald-type F-tests (s = 1) used instead.", call. = FALSE)
  }
  if(waldF == FALSE){
    thetaBetween <- sapply(H$HB, function(x){
      max(Re(eigen(solve(E$EB) %*% x)$values))
    })
    thetaWithin <- sapply(H$H, function(x){
      max(Re(eigen(solve(E$E) %*% x)$values))
    })

    pBetween <- sapply(H$HB, function(x){
      pBetween <- Matrix::rankMatrix(x + E$EB)
      attributes(pBetween) <- NULL
      pBetween
    })
    pWithin <- sapply(H$H, function(x){
      pWithin <- Matrix::rankMatrix(x + E$E)
      attributes(pWithin) <- NULL
      pWithin
    })

    X <- H$X
    rankX <- Matrix::rankMatrix(X)
    attributes(rankX) <- NULL
    v <- n - rankX
    Lbetween <- H$Lbetween
    L <- Lwithin <- NULL
  }else{
    HB <- H <- NULL
    X <- tt$Xbetween
    rankX <- Matrix::rankMatrix(X)
    attributes(rankX) <- NULL
    v <- n - rankX
    fMod <- mmemmuris::MLtoREML(fMod)

    L <- mmemmuris:::Lmatrix.ulm(fMod)
    Lwithin <- L[tt$within]
    if(length(tt$between) == 0L){
      L <- waldBetween <- X2Between <- thetaBetween <- NULL
    }else{
      term <- attr(terms(fMod), "term.labels")
      index <- grepl(paste0(":", tt$within[1], "|", tt$within[1], ":"), term)
      index <- unname(as.list(data.frame(rbind(tt$between, term[index]))))
      L <- lapply(index, function(x){ do.call("rbind", mmemmuris:::Lmatrix.ulm(fMod)[x]) })
      names(L) <- tt$between
      waldBetween <- mmemmuris::waldF(fMod, L)

      X2Between <- waldBetween$X2
      thetaBetween <- X2Between / v
    }

    waldWithin <- mmemmuris::waldF(fMod, Lwithin)
    X2Within <- waldWithin$X2
    thetaWithin <- X2Within / v

    Lbetween <- t(as.matrix(Matrix::expand(Matrix::lu(t(X) %*% X))$L))
    rownames(Lbetween) <- colnames(Lbetween) <- colnames(X)

    termsIndex <- attr(X, "assign")
    tableTerms <- table(termsIndex)
    termsLabels <- rep(c("(Intercept)", attr(terms(fMod), "term.labels")[attr(terms(fMod), "term.labels") %in% tt$between]), tableTerms)
    term <- c("(Intercept)", attr(terms(fMod), "term.labels")[attr(terms(fMod), "term.labels") %in% tt$between])
    rowIndex <- list()
    for(i in 1:length(term)){
      rowIndex[[i]] <- which(termsLabels == term[i])
      names(rowIndex)[i] <- term[i]
    }
    Lbetween <- lapply(rowIndex, function(x){
      Lbetween[x,, drop = FALSE]
    })
  }
  Qbetween <- Qwithin <- sapply(Lbetween, function(x){
    Q <- Matrix::rankMatrix(x %*% MASS::ginv(t(X) %*% X) %*% t(x))
    attributes(Q) <- NULL
    Q
  })

  if(waldF == FALSE){
    sBetween <- apply(cbind(pBetween, Qbetween), 1, min)
    Qwithin <- Qwithin[1:length(tt$within)]
  }else{
    if(length(tt$between) == 0L){
      Qbetween <- pBetween <- sBetween <- NULL
    }else{
      pqBetween <- waldBetween$NumDF
      Qbetween <- Qbetween[-1]
      pBetween <- pqBetween / Qbetween
      #names(pBetween) <- tt$between
      sBetween <- apply(cbind(pBetween, Qbetween), 1, min)
    }
    pqWithin <- waldWithin$NumDF
    Qwithin <- Qwithin[1:length(tt$within)]
    pWithin <- pqWithin / Qwithin
    names(pWithin) <- tt$within
  }
  sWithin <- apply(cbind(pWithin, Qwithin), 1, min)[1:length(tt$within)]

  theta <- list(EB = E$EB, HB = H$HB, thetaBetween = thetaBetween,
                sBetween = sBetween, pBetween = pBetween, qBetween = Qbetween,
                E = E$E, H = H$H,
                thetaWithin = thetaWithin, sWithin = sWithin,
                pWithin = pWithin, qWithin = Qwithin, X = X, v = v,
                Lbetween = Lbetween, Lwithin = Lwithin, L = L, n = n, type = type, waldF = waldF)

  class(theta) <- c("list", "theta.ulm")
  return(theta)
}

makeRoyGR.ulm <- function(theta){

  rWithin <- apply(cbind(theta$pWithin, theta$qWithin), 1, max)
  parmsRGRWithin <- data.frame(p = theta$pWithin, q = theta$qWithin, s = theta$sWithin,
                               r = rWithin, v = theta$v)
  FrgrWithin <- theta$thetaWithin * ((theta$v - rWithin + theta$qWithin) / rWithin)
  termsRGRWithin <- data.frame(dfNum = rWithin,
                               dfDen = theta$v - rWithin + theta$qWithin,
                               theta = theta$thetaWithin, F = FrgrWithin,
                               p.value = 1 - pf(FrgrWithin, rWithin, theta$v - rWithin + theta$qWithin))

  colnames(termsRGRWithin) <- c("Num df", "Den df", "theta", "F", "Pr(>F)")
  termsRGRWithin <- cbind(termsRGRWithin, mmemmuris::sigStars(termsRGRWithin$`Pr(>F)`))
  colnames(termsRGRWithin)[length(colnames(termsRGRWithin))] <- ""

  if(theta$waldF == TRUE){
    withinTests <- list(rgr = termsRGRWithin[theta$sWithin == 1, ],
                        parmsRGR = parmsRGRWithin,
                        approximation = "Pillai",
                        E = theta$E, H = theta$H,
                        Lwithin = theta$Lwithin, type = theta$type)
  }else{
    withinTests <- list(rgr = termsRGRWithin,
                        parmsRGR = parmsRGRWithin,
                        approximation = "Pillai",
                        E = theta$E, H = theta$H,
                        Lwithin = theta$Lwithin, type = theta$type)
  }

  if(!is.null(theta$thetaBetween)){
    rBetween <- apply(cbind(theta$pBetween, theta$qBetween), 1, max)
    parmsRGRBetween <- data.frame(p = theta$pBetween, q = theta$qBetween,
                                  s = theta$sBetween, r = rBetween,
                                  v = theta$v)
    FrgrBetween <- theta$thetaBetween * ((theta$v - rBetween + theta$qBetween) / rBetween)
    termsRGRBetween <- data.frame(dfNum = rBetween,
                                  dfDen = theta$v - rBetween + theta$qBetween,
                                  theta = theta$thetaBetween, F = FrgrBetween,
                                  p.value = 1 - pf(FrgrBetween, rBetween, theta$v - rBetween + theta$qBetween))

    colnames(termsRGRBetween) <- c("Num df", "Den df", "theta", "F", "Pr(>F)")
    termsRGRBetween <- cbind(termsRGRBetween, mmemmuris::sigStars(termsRGRBetween$`Pr(>F)`))
    colnames(termsRGRBetween)[length(colnames(termsRGRBetween))] <- ""

    if(theta$waldF){
      betweenTests <- list(rgr = termsRGRBetween[rownames(termsRGRBetween) != "(Intercept)" & theta$sBetween == 1, ],
                           parmsRGR = parmsRGRBetween[rownames(parmsRGRBetween) != "(Intercept)", ],
                           approximation = "Pillai",
                           EB = theta$EB, HB = theta$HB,
                           Lbetween = theta$Lbetween, L = theta$L, type = theta$type,
                           Vmatrix = theta$EB / theta$n)
    }else{
      betweenTests <- list(rgr = termsRGRBetween[rownames(termsRGRBetween) != "(Intercept)", ],
                           parmsRGR = parmsRGRBetween[rownames(parmsRGRBetween) != "(Intercept)", ],
                           approximation = "Pillai",
                           EB = theta$EB, HB = theta$HB,
                           Lbetween = theta$Lbetween, L = theta$L, type = theta$type,
                           Vmatrix = theta$EB / theta$n)
    }

    if(nrow(betweenTests$rgr) > 0L)
      rgr <- list(betweenTests = betweenTests, withinTests = withinTests,
                  waldF = theta$waldF)
    else
      rgr <- list(withinTests = withinTests,
                  waldF = theta$waldF)
    class(rgr) <- c("list", "multiTest", "rgr", "ulm")
  }else{
    rgr <- list(withinTests = withinTests,
                waldF = theta$waldF)
    class(rgr) <- c("list", "multiTest", "rgr", "ulm")
  }
  return(rgr)
}

#' Likelihood Ratio Test for Nonsymbolically Nested Models (Compound Symmetry Covariance Nested Within an Unstructured Covariance)
#'
#' @export
#'
#' @description
#' This function will provide a possible rationale for choosing to use the
#' multivariate tests (unstructured) or univariate tests (compound symmetry)
#' based on a likelihood ratio test, AIC, or BIC of the covariance structures
#' for a \code{\link[nlme]{gls}}, \code{\link[nlme]{lme}}, or \code{mlm} model
#' object.
#'
#' @details
#'
#' **Marginal and Mixed-Effects Models**
#'
#' The models with compound symmetry and unstructured covariance structures are
#' fit with maximum likelihood so that a likelihood ratio test can be done for
#' the covariance structures.  The maximum log-likelihood functions for the
#' models with compound symmetry and unstructured covariances, respectively, are
#'
#' \eqn{logLik(M_R) = (-1 / 2)log(det(}**V_{CS}** \eqn{)) - ((1 / 2)}**r'V_{CS}** \eqn{^ {-1}}**r** \eqn{) - (N / 2)log(2\pi)}
#'
#' and
#'
#' \eqn{logLik(M_F) = (-1 / 2)log(det(}**V_{UN}** \eqn{)) - ((1 / 2)}**r'V_{UN}** \eqn{^ {-1}}**r** \eqn{) - (N / 2)log(2\pi)}
#'
#' where **V** is the marginal variance-covariance matrix of the random/repeated
#' effects, **r** is the vector of model residuals, and \eqn{N} is the number of
#' observations of data in the long format.
#'
#' The likelihood ratio (LR) statistic approximately follows a \eqn{\chi ^ 2}
#' distribution \deqn{LR = -2logLik(M_R) - {-2logLik(M_F)}} with
#' \eqn{df_{M_R} - df_{M_F}} degrees of freedom.
#'
#' The Akaike's Information Criterion (AIC) is defined as
#' \deqn{AIC = 2df_M - 2logLik(M)}
#' where \eqn{df_M} is the number of parameter estimates in the model.
#'
#' The Bayesian Information Criterion (BIC) is defined as
#' \deqn{BIC = df_Mlog(N) - 2logLik(M)}
#' where \eqn{df_M} is the number of parameter estimates in the model.
#'
#' **MANOVA/RM ANOVA**
#'
#' For MANOVA/RM ANOVA, the (log-)likelihoods are not computed since method of
#' moments is used instead of (RE)ML.  Therefore, this could be viewed as a fake
#' solution used for illustration purposes.  We start by obtaining the marginal
#' variance-covariance matrix of the the random/repeated effects from the
#' between-subject effects **E** matrix.  The model with an unstructured
#' covariance using method of moments coincides with the maximum likelihood
#' estimates for the marginal variance-covariance matrix.  However, the model
#' with a compound symmetry covariance using method of moments coincides with
#' the restricted maximum likelihood estimates for the marginal
#' variance-covariance matrix.  Since the maximum likelihood estimates are
#' biased, we need to undo Bessel's correction by multiplying the marginal
#' variance-covariance matrix "obtained" by REML by \eqn{v / n} where \eqn{v} is
#' the residual degrees of freedom for the multivariate linear model and \eqn{n}
#' is the total number of observations in the dataset in the wide format.
#'
#' The maximum log-likelihood functions for the
#' models with compound symmetry and unstructured covariances, respectively, are
#'
#' \eqn{logLik(M_R) = (-1 / 2)log(det(}**V_{CS}** \eqn{)) - ((1 / 2)}**r'V_{CS}** \eqn{^ {-1}}**r** \eqn{) - (N / 2)log(2\pi)}
#'
#' and
#'
#' \eqn{logLik(M_F) = (-1 / 2)log(det(}**V_{UN}** \eqn{)) - ((1 / 2)}**r'V_{UN}** \eqn{^ {-1}}**r** \eqn{) - (N / 2)log(2\pi)}
#'
#' where **V** is the marginal variance-covariance matrix of the random/repeated
#' effects, **r** is the vector of model residuals, and \eqn{N} is the number of
#' observations of data in the long format.
#'
#' The likelihood ratio (LR) statistic approximately follows a \eqn{\chi ^ 2}
#' distribution \deqn{LR = -2logLik(M_R) - {-2logLik(M_F)}} with
#' \eqn{df_{M_R} - df_{M_F}} degrees of freedom.
#'
#' The Akaike's Information Criterion (AIC) is defined as
#' \deqn{AIC = 2df_M - 2logLik(M)}
#' where \eqn{df_M} is the number of parameter estimates in the model.
#'
#' The Bayesian Information Criterion (BIC) is defined as
#' \deqn{BIC = df_Mlog(N) - 2logLik(M)}
#' where \eqn{df_M} is the number of parameter estimates in the model.
#'
#' @param fMod A model of class \code{\link[nlme]{gls}}, \code{\link[nlme]{lme}},
#' or \code{mlm}.
#'
#' @references {\url{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.3/statug/statug_mixed_details01.htm#statug.mixed.mixedmmtgandr}}
#'
#' @returns \tabular{ll}{
#'    `df` \tab The degrees of freedom used for the models. \cr
#'    \tab \cr
#'    `AIC` \tab The Akaike information criterions (AICs) for the models. \cr
#'    \tab \cr
#'    `BIC` \tab The Bayesian information criterions (BICs) for the models. \cr
#'    \tab \cr
#'    `logLik` \tab The log-likelihoods for the models. \cr
#'    \tab \cr
#'    `L.Ratio` \tab The likelihood ratio for the models. \cr
#'    \tab \cr
#'    `Pr(>Chisq)` \tab The p-value for the likelihood ratio test for the models. \cr
#' }
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#' nestedCovariance(glsMod)
#'
#' # Linear mixed-effects model
#' # Unstructured covariance
#' lmeMod <- lme(Distance ~ Sex * Age,
#'               random = ~ Age | Subject,
#'               na.action = na.omit,
#'               data = orthodontLong)
#' nestedCovariance(lmeMod)
#'
#' # Multivariate linear model
#' manovaMod <- manova(cbind(Distance8, Distance10, Distance12, Distance14) ~ Sex,
#'                     data = orthodontWide)
#' nestedCovariance(manovaMod)
#'
#' @seealso
#' {\code{\link[mmemmuris]{covStruct}}, \code{\link[mmemmuris]{REMLtoML}}}

nestedCovariance <- function(fMod){
  if(any(class(fMod) %in% c("gls", "lme"))){
    if(mmemmuris::covStruct(fMod) %in% c("cs", "other"))
      stop("Only unstructured covariance models are allowed.", call. = FALSE)
    # tt <- mmemmuris::termsType(fMod)
    # data <- mmemmuris:::getDataset(fMod)
    # Refit models w/ ML if fit by REML
    fMod <- mmemmuris::REMLtoML(fMod)
    if(class(fMod) == "gls"){
      rMod <- update(fMod,
                     correlation = nlme::corCompSymm(form = as.formula(paste("~ 1 | ", as.character(nlme::getGroupsFormula(fMod))[2]))),
                     weights = NULL,
                     na.action = na.omit)
    }else{
      # Update the model to compound symmetry covariance
      rMod <- update(fMod,
                     random = as.formula(paste("~ 1 | ", as.character(nlme::getGroupsFormula(fMod))[2], "")),
                     correlation = NULL,
                     weights = NULL,
                     na.action = na.omit)
    }
    nestedCov <- anova(rMod, fMod, test = TRUE)
    class(nestedCov) <- "data.frame"
    rownames(nestedCov)[2] <- "Unstructured"
    if(class(fMod) %in% "gls")
      rownames(nestedCov)[1] <- "Compound Symmetry"
    else
      rownames(nestedCov)[1] <- "Random Intercepts"
    nestedCov <- nestedCov[, !colnames(nestedCov) %in% c("call", "Model", "Test")]
  }else if(any(class(fMod) %in% "mlm")){
    dat <- na.omit(mmemmuris:::getDataset(fMod))
    V <- mmemmuris::Vmatrix(fMod)$unV
    V <- as.matrix(bdiag(replicate(nrow(dat), V, simplify = FALSE)))
    r <- cbind(c(t(residuals(fMod))))
    N <- nrow(dat) * length(colnames(coef(fMod)))
    UN <- (-1 / 2) * log(det(V)) - ((1 / 2) * t(r) %*% solve(V) %*% r) - (N / 2) * log(2 * pi)
    V <- (fMod$df.residual / nrow(dat)) * mmemmuris::Vmatrix(fMod)$csV
    V <- as.matrix(bdiag(replicate(nrow(dat), V, simplify = FALSE)))
    CS <- (-1 / 2) * log(det(V)) - ((1 / 2) * t(r) %*% solve(V) %*% r) - (N / 2) * log(2 * pi)
    dfModel <- c(length(coef(fMod)) + 2, length(coef(fMod)) + ((length(colnames(coef(fMod))) * (length(colnames(coef(fMod))) + 1)) / 2))
    X2 <- -2 * CS - (-2 * UN)
    lik <- c(CS, UN)
    nestedCov <- data.frame(df = dfModel,
                            AIC = 2 * dfModel - 2 * lik,
                            BIC = dfModel * log(N) - 2 * lik,
                            logLik = lik,
                            L.Ratio = c(NA, X2),
                            p.value = c(NA, 1 - pchisq(X2, diff(dfModel))))
    rownames(nestedCov) <- c("Compound Symmetry", "Unstructured")
  }else{
    stop('Please input a model of class "gls", "lme", or "mlm".', call. = FALSE)
  }
  colnames(nestedCov)[6] <- "Pr(>Chisq)"
  nestedCov <- cbind(nestedCov, mmemmuris::sigStars(nestedCov$`Pr(>Chisq)`))
  nestedCov <- list(nestedCov = nestedCov)
  class(nestedCov) <- c("list", "nestedCovariance")
  return(nestedCov)
}

#' Residual Denominator Degrees of Freedom
#'
#' @export
#'
#' @description
#' This function calculates the residual denominator degrees of freedom for marginal
#' (\code{\link[nlme]{gls}}) and mixed-effects (\code{\link[nlme]{lme}}) models.
#'
#' @param fMod A model of class \code{\link[nlme]{gls}} or \code{\link[nlme]{lme}}.
#'
#' @details
#' The residual denominator degrees of freedom are calculated as
#'
#' \eqn{ddf = N - rank(}**X**\eqn{)}
#'
#' where \eqn{N} is the total number of observations of data in the long format
#' and **X** is the design matrix of the independent variables.
#'
#' @references {\url{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.3/statug/statug_glimmix_details38.htm}}
#'
#' @returns
#' \tabular{ll}{
#'    `between`  \tab A character vector containing the between-subject effects
#'    for the model. \cr
#'    `within`  \tab A character vector containing the within-subject effects
#'    for the model. \cr
#'    ``Den df``  \tab A vector containing the denominator degrees of freedom
#'    for each term in the model. \cr
#'    `Xbetween`  \tab The between-subject effects design matrix. \cr
#'    `X`  \tab The fixed effects design matrix. \cr
#' \tab \cr
#' }
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#' ddfResidual(glsMod)
#' # Compound symmetry covariance
#' glsModCS <- gls(Distance ~ Sex * Age,
#'                 correlation = corCompSymm(form = ~ 1 | Subject),
#'                 na.action = na.omit,
#'                 data = orthodontLong)
#' ddfResidual(glsModCS)

ddfResidual <- function(fMod){
  if(any(class(fMod) %in% c("gls", "lme"))){
    tt <- mmemmuris::termsType(fMod)
    X <- model.matrix(fMod, data = mmemmuris:::getDataset(fMod))
    rDenDF <- fMod$dims$N - Matrix::rankMatrix(X)
    attributes(rDenDF) <- NULL
    residualDenDF <- rep(rDenDF, length(c(tt$between, tt$within)))
    names(residualDenDF) <- c(tt$between, tt$within)
    return(list(between = tt$between, within = tt$within,
                `Den df` = residualDenDF,
                Xbetween = tt$Xbetween, X = X))
  }else
    stop('The model must be a "gls" or "lme" object.', call. = FALSE)
}

#' Between-Within Denominator Degrees of Freedom
#'
#' @export
#'
#' @description
#' This function calculates the between-within denominator degrees of freedom for
#' marginal (\code{\link[nlme]{gls}}) and mixed-effects (\code{\link[nlme]{lme}}) models.
#'
#' @param fMod A model of class \code{\link[nlme]{gls}} or \code{\link[nlme]{lme}}.
#'
#' @details
#' The residual denominator degrees of freedom are calculated as
#'
#' \eqn{ddf = N - rank(}**X**\eqn{)}
#'
#' where \eqn{N} is the total number of observations of data in the long format
#' and **X** is the design matrix of the independent variables.
#'
#' The denominator degrees of freedom are equal to
#'
#' \eqn{ddfBetween = n - rank(}**X_b** \eqn{)}
#'
#' where \eqn{n} is the total number of observations of data in the wide format
#' and **X_b** is the design matrix of the between-subject independent
#' variables in the wide format
#'
#' for between-subject effects and
#'
#' \eqn{ddfWithin = ddf - ddfBetween}
#'
#' for within-subject effects.
#'
#' @references {\url{https://documentation.sas.com/doc/en/statug/15.2/statug_glimmix_details37.htm}}
#'
#' @returns
#' \tabular{ll}{
#'    `between`  \tab A character vector containing the between-subject effects
#'    for the model. \cr
#'    `within`  \tab A character vector containing the within-subject effects
#'    for the model. \cr
#'    `bDenDF`  \tab A vector containing the between-subject degrees of freedom. \cr
#'    `wDenDF`  \tab A vector containing the within-subject degrees of freedom. \cr
#'    `rDenDF`  \tab A vector containing the residual degrees of freedom. \cr
#'    ``Den df``  \tab A vector containing the denominator degrees of freedom
#'    for each term in the model. \cr
#'    `Xbetween`  \tab The between-subject effects design matrix. \cr
#' \tab \cr
#' }
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#' ddfBW(glsMod)
#' # Compound symmetry covariance
#' glsModCS <- gls(Distance ~ Sex * Age,
#'                 correlation = corCompSymm(form = ~ 1 | Subject),
#'                 na.action = na.omit,
#'                 data = orthodontLong)
#' ddfBW(glsModCS)
#'
#' @seealso
#' {\code{\link[mmemmuris]{ddfResidual}}}

ddfBW <- function(fMod){
  if(any(class(fMod) %in% c("gls", "lme"))){
    residualDenDF <- mmemmuris::ddfResidual(fMod)
    if(class(fMod) %in% "gls"){
      bDenDF <- length(levels(getGroups(fMod))) - Matrix::rankMatrix(residualDenDF$Xbetween)
      attributes(bDenDF) <- NULL
      if("corSymm" %in% class(fMod$modelStruct$corStruct)){
        bwDenDF <- rep(bDenDF, length(c(residualDenDF$between, residualDenDF$within)))
      }else if("corCompSymm" %in% class(fMod$modelStruct$corStruct)){
        wDenDF <- residualDenDF$`Den df`[residualDenDF$within] - bDenDF
        bwDenDF <- c(rep(bDenDF, length(residualDenDF$between)), wDenDF)
      }
    }else if(class(fMod) %in% "lme"){
      bDenDF <- length(levels(getGroups(fMod))) - Matrix::rankMatrix(residualDenDF$Xbetween)
      attributes(bDenDF) <- NULL
      wDenDF <- residualDenDF$`Den df`[residualDenDF$within] - bDenDF
      bwDenDF <- c(rep(bDenDF, length(residualDenDF$between)), wDenDF)
    }
    names(bwDenDF) <- c(residualDenDF$between, residualDenDF$within)
    if(class(fMod) %in% "gls"){
      if(exists("wDenDF"))
        return(list(between = residualDenDF$between, within = residualDenDF$within,
                    bDenDF = bDenDF, wDenDF = unname(wDenDF[1]),
                    rDenDF = unname(residualDenDF$`Den df`[1]),
                    `Den df` = bwDenDF,
                    Xbetween = residualDenDF$Xbetween))
      else
        return(list(between = residualDenDF$between, within = residualDenDF$within,
                    bDenDF = bDenDF,
                    rDenDF = unname(residualDenDF$`Den df`[1]),
                    `Den df` = bwDenDF,
                    Xbetween = residualDenDF$Xbetween))
    }else{
      return(list(between = residualDenDF$between, within = residualDenDF$within,
                  bDenDF = bDenDF, wDenDF = unname(wDenDF[1]),
                  rDenDF = unname(residualDenDF$`Den df`[1]),
                  `Den df` = bwDenDF,
                  Xbetween = residualDenDF$Xbetween))
    }
  }else
    stop('Please provide a model of class "gls" or "lme".', call. = FALSE)
}

#' Significance Level Stars for p-values
#'
#' @export
#' @description
#' This creates a vector of significance level stars for p-values.
#'
#' Significance level codes: \tabular{lll}{ \cr
#' \[0, 0.001\]    \tab - \tab *** \cr
#' (0.001, 0.01] \tab - \tab ** \cr
#' (0.01, 0.05]  \tab - \tab * \cr
#' (0.05, 0.1]   \tab - \tab . \cr
#' (0.1, 1]      \tab -
#' }
#'
#' @param pValues A vector of p-values.
#'
#' @returns
#' This function returns a vector of significance stars.
#'
#' @examples
#' pValues <- c(0, 1e-16, 0.005, 0.025, 0.075, 0.5, 1)
#' sigStars(pValues)

sigStars <- function(pValues){
  if(!is.numeric(na.omit(pValues)) | any(na.omit(pValues) < 0 | na.omit(pValues) > 1))
    stop("p-values must be numeric and bounded in [0, 1].", call. = FALSE)
  sig <- ifelse(0 <= pValues & pValues <= 0.001, "***",
                ifelse(0.001 < pValues & pValues <= 0.01, "**",
                       ifelse(0.01 < pValues & pValues <= 0.05, "*",
                              ifelse(0.05 < pValues & pValues <= 0.1, ".",
                                     ifelse(0.1 < pValues & pValues <= 1, " ", NA)))))
  return(sig)
}

#' Print a Univariate Tests Object
#'
#' @exportS3Method
print.uniTest <- function(uniTest, V = FALSE, H = FALSE, E = FALSE,
                          digits = 4){
  if(digits < 3)
    stop("No less than 3 digits can be printed.", call. = FALSE)
  cat("\n----------------\n")
  cat("Univariate Tests")
  cat("\n----------------\n")
  if(V){
    cat("\nMarginal Variance-Covariance Matrix Assuming Sphericity\n")
    cat("-------------------------------------------------------\n")
    Vmatrix <- as.data.frame(format(round(uniTest$csV$V, digits), nsmall = digits,
                                    scientific = FALSE))
    colnames(Vmatrix) <- NULL
    print(Vmatrix, row.names = FALSE)
    cat("\n")
    cat("\nMarginal Correlation Matrix Assuming Sphericity\n")
    cat("-----------------------------------------------\n")
    Vcorr <- as.data.frame(format(round(cov2cor(uniTest$csV$V), digits), nsmall = digits,
                                  scientific = FALSE))
    colnames(Vcorr) <- NULL
    print(Vcorr, row.names = FALSE)
    cat("\n")
  }
  if(nrow(uniTest$betweenTests$betweenTests[rownames(uniTest$betweenTests$betweenTests) != "Residuals  ", ]) > 0){
    cat("\n")
    type <- paste0("Type ", uniTest$betweenTests$type, " Tests: Between-Subject Effects")
    cat(type)
    cat("\n")
    cat(paste(rep("-", nchar(type)), collapse = ""), "\n")
    betweenTests <- cbind(format(round(uniTest$betweenTests$betweenTests[, -ncol(uniTest$betweenTests$betweenTests)], digits), nsmall = digits,
                                 scientific = FALSE),
                          uniTest$betweenTests$betweenTests[, ncol(uniTest$betweenTests$betweenTests)])
    colnames(betweenTests)[length(colnames(betweenTests))] <- ""
    print(betweenTests)
  }
  if(H & ("mlm" %in% class(uniTest))){
    cat("\n")
    type <- paste0("Type ", uniTest$withinTests$type, ": Within-Subject SSCP H matrix")
    cat(type)
    cat("\n")
    cat(paste(rep("-", nchar(type)), collapse = ""), "\n")
    print(lapply(uniTest$H, function(x){
      H <- as.data.frame(format(round(x, digits), nsmall = digits,
                                scientific = FALSE))
      colnames(H) <- NULL
      return(H)
    }), row.names = FALSE)
  }
  if(E){
    cat("\nWithin-Subject SSCP E matrix\n")
    cat("----------------------------\n")
    E <- as.data.frame(format(round(uniTest$E, digits), nsmall = digits,
                              scientific = FALSE))
    colnames(E) <- NULL
    print(E, row.names = FALSE)
  }
  if(!is.null(uniTest$sphericityTests))
    print(uniTest$sphericityTests, digits = digits)
  else{
    # cat("\n")
    # warning("The SSCP E matrix is singular or nearly singular or sphericity is automatically satisfied.  As a result, sphericity tests and epsilon corrections are not printed.  Try adjusting the tolerance at your own risk if desired.", call. = FALSE)
  }
  print(uniTest$withinTests, digits = digits)
  if(!is.null(uniTest$epsCorrect))
    print(uniTest$epsCorrect, digits = digits)
}

#' Print a Wilks' Lambda Object
#'
#' @exportS3Method
print.wilks <- function(wilks, digits = 4){
  # if(wilks$singular == FALSE){
    if(digits < 3)
      stop("No less than 3 digits can be printed.", call. = FALSE)
    if(length(wilks$betweenTests$wilks) > 0L){
      if(nrow(wilks$betweenTests$wilks) > 0L){
        cat("\n")
        type <- paste0("Type ", wilks$betweenTests$type, " Tests: Wilks' Lambda for Between-Subject Effects")
        cat(type)
        cat("\n")
        cat(paste(rep("-", nchar(type)), collapse = ""), "\n")
        print(wilks$betweenTests$parmsWilks)
        cat("\n")
        wilksBetween <- cbind(format(round(wilks$betweenTests$wilks[, -ncol(wilks$betweenTests$wilks)], digits), nsmall = digits,
                                     scientific = FALSE),
                              wilks$betweenTests$wilks[, ncol(wilks$betweenTests$wilks)])
        colnames(wilksBetween)[length(colnames(wilksBetween))] <- ""
        print(wilksBetween)
        cat("---")
        cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
        cat(paste0("Approximation method: ", wilks$betweenTests$approximation, "\n\n"))
        if(wilks$betweenTests$approximation == "Rao"){
          for(i in 1:nrow(wilks$betweenTests$parmsWilks)){
            if(wilks$betweenTests$parmsWilks$s[i] <= 2L)
              cat("The F-statistic for", rownames(wilks$betweenTests$parmsWilks)[i], "is exact.\n")
          }
        }
      }
    }
    cat("\n")
    type <- paste0("Type ", wilks$withinTests$type, " Tests: Wilks' Lambda for Within-Subject Effects")
    cat(type)
    cat("\n")
    cat(paste(rep("-", nchar(type)), collapse = ""), "\n")
    print(wilks$withinTests$parmsWilks)
    cat("\n")
    wilksWithin <- cbind(format(round(wilks$withinTests$wilks[, -ncol(wilks$withinTests$wilks)], digits), nsmall = digits,
                                scientific = FALSE),
                         wilks$withinTests$wilks[, ncol(wilks$withinTests$wilks)])
    colnames(wilksWithin)[length(colnames(wilksWithin))] <- ""
    print(wilksWithin)
    cat("---")
    cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
    cat(paste0("Approximation method: ", wilks$withinTests$approximation, "\n\n"))
    if(wilks$withinTests$approximation == "Rao"){
      for(i in 1:nrow(wilks$withinTests$parmsWilks)){
        if(wilks$withinTests$parmsWilks$s[i] <= 2)
          cat("The F-statistic for", rownames(wilks$withinTests$parmsWilks)[i], "is exact.\n")
      }
    }
    cat("\n")
}

#' Print a Pillai's Trace Object
#'
#' @exportS3Method
print.pillai <- function(pillai, digits = 4){
  # if(pillai$singular == FALSE){
  if(digits < 3)
    stop("No less than 3 digits can be printed.", call. = FALSE)
  if(!is.null(pillai$betweenTests)){
    cat("\n")
    type <- paste0("Type ", pillai$betweenTests$type, " Tests: Pillai's Trace for Between-Subject Effects")
    cat(type)
    cat("\n")
    cat(paste(rep("-", nchar(type)), collapse = ""), "\n")
    print(pillai$betweenTests$parmsPillai)
    cat("\n")
    pillaiBetween <- cbind(format(round(pillai$betweenTests$pillai[, -ncol(pillai$betweenTests$pillai)], digits), nsmall = digits,
                                  scientific = FALSE),
                           pillai$betweenTests$pillai[, ncol(pillai$betweenTests$pillai)])
    colnames(pillaiBetween)[length(colnames(pillaiBetween))] <- ""
    print(pillaiBetween)
    cat("---")
    cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
    cat(paste0("Approximation method: ", pillai$betweenTests$approximation, "\n\n"))
    if(nrow(pillai$betweenTests$parmsPillai) > 0L){
      for(i in 1:nrow(pillai$betweenTests$parmsPillai)){
        if(pillai$betweenTests$parmsPillai$s[i] == 1L)
          cat("The F-statistic for", rownames(pillai$betweenTests$parmsPillai)[i], "is exact.\n")
      }
    }
  }
  cat("\n")
  type <- paste0("Type ", pillai$withinTests$type, " Tests: Pillai's Trace for Within-Subject Effects")
  cat(type)
  cat("\n")
  cat(paste(rep("-", nchar(type)), collapse = ""), "\n")
  print(pillai$withinTests$parmsPillai)
  cat("\n")
  pillaiWithin <- cbind(format(round(pillai$withinTests$pillai[, -ncol(pillai$withinTests$pillai)], digits), nsmall = digits,
                               scientific = FALSE),
                        pillai$withinTests$pillai[, ncol(pillai$withinTests$pillai)])
  colnames(pillaiWithin)[length(colnames(pillaiWithin))] <- ""
  print(pillaiWithin)
  cat("---")
  cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat(paste0("Approximation method: ", pillai$withinTests$approximation, "\n\n"))
  if(nrow(pillai$withinTests$parmsPillai) > 0L){
    for(i in 1:nrow(pillai$withinTests$parmsPillai)){
      if(pillai$withinTests$parmsPillai$s[i] == 1L)
        cat("The F-statistic for", rownames(pillai$withinTests$parmsPillai)[i], "is exact.\n")
    }
  }
  cat("\n")
}

#' Print a Hotelling-Lawley Trace Object
#'
#' @exportS3Method
print.hlt <- function(hlt, digits = 4){
  if(digits < 3)
    stop("No less than 3 digits can be printed.", call. = FALSE)
  if(length(hlt$betweenTests$hlt) > 0L){
    if(nrow(hlt$betweenTests$hlt) > 0L){
      cat("\n")
      type <- paste0("Type ", hlt$betweenTests$type, " Tests: Hotelling-Lawley Trace for Between-Subject Effects")
      cat(type)
      cat("\n")
      cat(paste(rep("-", nchar(type)), collapse = ""), "\n")
      print(hlt$betweenTests$parmsHL)
      cat("\n")
      hltBetween <- cbind(format(round(hlt$betweenTests$hlt[, -ncol(hlt$betweenTests$hlt)], digits), nsmall = digits,
                                 scientific = FALSE),
                          hlt$betweenTests$hlt[, ncol(hlt$betweenTests$hlt)])
      colnames(hltBetween)[length(colnames(hltBetween))] <- ""
      print(hltBetween)
      cat("---")
      cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
      cat(paste0("Approximation method: ", hlt$betweenTests$approximation, "\n\n"))
      if(hlt$betweenTests$approximation %in% c("McKeon", "Pillai-Samson")){
        for(i in 1:nrow(hlt$betweenTests$parmsHL)){
          if(hlt$betweenTests$parmsHL$s[i] == 1)
            cat("The F-statistic for", rownames(hlt$betweenTests$parmsHL)[i], "is exact.\n")
        }
      }
    }
  }
  cat("\n")
  type <- paste0("Type ", hlt$withinTests$type, " Tests: Hotelling-Lawley Trace for Within-Subject Effects")
  cat(type)
  cat("\n")
  cat(paste(rep("-", nchar(type)), collapse = ""), "\n")
  print(hlt$withinTests$parmsHL)
  cat("\n")
  hltWithin <- cbind(format(round(hlt$withinTests$hlt[, -ncol(hlt$withinTests$hlt)], digits), nsmall = digits,
                            scientific = FALSE),
                     hlt$withinTests$hlt[, ncol(hlt$withinTests$hlt)])
  colnames(hltWithin)[length(colnames(hltWithin))] <- ""
  print(hltWithin)
  cat("---")
  cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat(paste0("Approximation method: ", hlt$withinTests$approximation, "\n\n"))
  if(hlt$withinTests$approximation %in% c("McKeon", "Pillai-Samson")){
    for(i in 1:nrow(hlt$withinTests$parmsHL)){
      if(hlt$withinTests$parmsHL$s[i] == 1)
        cat("The F-statistic for", rownames(hlt$withinTests$parmsHL)[i], "is exact.\n")
    }
  }
  cat("\n")
}

#' Print a Roy's Greatest Root Object
#'
#' @exportS3Method
print.rgr <- function(rgr, digits = 4){
  if(digits < 3)
    stop("No less than 3 digits can be printed.", call. = FALSE)
  if(!is.null(rgr$betweenTests$rgr)){
    cat("\n")
    type <- paste0("Type ", rgr$betweenTests$type, " Tests: Roy's Greatest Root for Between-Subject Effects")
    cat(type)
    cat("\n")
    cat(paste(rep("-", nchar(type)), collapse = ""), "\n")
    print(rgr$betweenTests$parmsRGR)
    cat("\n")
    rgrBetween <- cbind(format(round(rgr$betweenTests$rgr[, -ncol(rgr$betweenTests$rgr)], digits), nsmall = digits,
                               scientific = FALSE),
                        rgr$betweenTests$rgr[, ncol(rgr$betweenTests$rgr)])
    colnames(rgrBetween)[length(colnames(rgrBetween))] <- ""
    print(rgrBetween)
    cat("---")
    cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
    cat(paste0("Approximation method: ", rgr$betweenTests$approximation, "\n\n"))
    if(nrow(rgr$betweenTests$parmsRGR) > 0L){
      for(i in 1:nrow(rgr$betweenTests$parmsRGR)){
        if(rgr$betweenTests$parmsRGR$s[i] <= 1)
          cat("The F-statistic for", rownames(rgr$betweenTests$parmsRGR)[i], "is exact.\n")
        else if(rgr$waldF == FALSE)
          cat("The F-statistic for", rownames(rgr$betweenTests$parmsRGR)[i], "is an upper bound.\n")
      }
    }
  }
  cat("\n")
  type <- paste0("Type ", rgr$withinTests$type, " Tests: Roy's Greatest Root for Within-Subject Effects")
  cat(type)
  cat("\n")
  cat(paste(rep("-", nchar(type)), collapse = ""), "\n")
  print(rgr$withinTests$parmsRGR)
  cat("\n")
  rgrWithin <- cbind(format(round(rgr$withinTests$rgr[, -ncol(rgr$withinTests$rgr)], digits), nsmall = digits,
                            scientific = FALSE),
                     rgr$withinTests$rgr[, ncol(rgr$withinTests$rgr)])
  colnames(rgrWithin)[length(colnames(rgrWithin))] <- ""
  print(rgrWithin)
  cat("---")
  cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat(paste0("Approximation method: ", rgr$withinTests$approximation, "\n\n"))
  if(nrow(rgr$withinTests$parmsRGR) > 0L){
    for(i in 1:nrow(rgr$withinTests$parmsRGR)){
      if(rgr$withinTests$parmsRGR$s[i] <= 1)
        cat("The F-statistic for", rownames(rgr$withinTests$parmsRGR)[i], "is exact.\n")
      else if(rgr$waldF == FALSE)
        cat("The F-statistic for", rownames(rgr$withinTests$parmsRGR)[i], "is an upper bound.\n")
    }
  }
  cat("\n")
}

#' Print a Within Tests Object
#'
#' @exportS3Method
print.withinTests <- function(withinTests, digits = 4){
  if(digits < 3)
    stop("No less than 3 digits can be printed.", call. = FALSE)
  cat("\n")
  type <- paste0("Type ", withinTests$type, " Tests: Within-Subject Effects Assuming Sphericity")
  cat(type)
  cat("\n")
  cat(paste(rep("-", nchar(type)), collapse = ""), "\n")
  within <- cbind(format(round(withinTests$withinTests[, -ncol(withinTests$withinTests)], digits), nsmall = digits,
                          scientific = FALSE),
                  withinTests$withinTests[, ncol(withinTests$withinTests)])
  colnames(within)[length(colnames(within))] <- ""
  print(within)
  cat("---")
  cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat(paste0("Denominator degrees of freedom: ", withinTests$ddf, "\n\n"))
}

#' Print an Epsilon Corrections Object
#'
#' @exportS3Method
print.epsCorrect <- function(epsCorrect, digits = 4){
  if(digits < 3)
    stop("No less than 3 digits can be printed.", call. = FALSE)
  correctMethod <- epsCorrect$correctMethod
  epsCorrect <- cbind(format(round(epsCorrect$epsCorrect[, -ncol(epsCorrect$epsCorrect)],
                                   digits), nsmall = digits, scientific = FALSE),
                      epsCorrect$epsCorrect[, ncol(epsCorrect$epsCorrect)])
  colnames(epsCorrect)[length(colnames(epsCorrect))] <- ""
  cat("\nEpsilon Correction for a Violation of Sphericity\n")
  cat("------------------------------------------------\n")
  print(epsCorrect)
  cat("---")
  cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  cat(paste0("Epsilon correction: ", correctMethod, "\n\n"))
}

#' Print a Nested Covariance Object
#'
#' @exportS3Method
print.nestedCovariance <- function(nestedCovariance, digits = 4){
  if(digits < 3)
    stop("No less than 3 digits can be printed.", call. = FALSE)
  cat("\nLikelihood Ratio Test for Nested Covariance Structures\n")
  cat("------------------------------------------------------\n")
  nestedCov <- cbind(format(round(nestedCovariance$nestedCov[, -ncol(nestedCovariance$nestedCov)],
                                     digits), nsmall = digits, scientific = FALSE),
                     nestedCovariance$nestedCov[, ncol(nestedCovariance$nestedCov)])
  colnames(nestedCov)[length(colnames(nestedCov))] <- nestedCov[1, length(colnames(nestedCov))] <- ""
  print(nestedCov)
  cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}

#' Create the L coefficient contrast matrix
#'
#' @description
#' This function creates a Type 1 Sum of Squares **L** matrix.
#'
#' @details
#' The contrast coefficients matrix (**L**) is created from an \eqn{LU}
#' decomposition of the crossproduct matrix **X'X** with the design matrix **X**
#' such that
#' \deqn{PA = LU}
#'
#' The **L** matrix is then defined as
#'
#' **L** \eqn{= L'}

Lmatrix.ulm <- function(fMod){
  terms <- terms(fMod)
  tt <- mmemmuris::termsType(fMod)
  X <- model.matrix(fMod, data = mmemmuris:::getDataset(fMod))
  if (ncol(X) == 1L)
    matrix(1L)
  else
    L <- t(as.matrix(Matrix::expand(Matrix::lu(t(X) %*% X))$L))
  rownames(L) <- colnames(L) <- colnames(X)
  termsIndex <- attr(X, "assign")
  tableTerms <- table(termsIndex)
  if(names(tableTerms)[1] == "0")
    termsLabels <- rep(c("(Intercept)", attr(terms, "term.labels")), tableTerms)
  else
    termsLabels <- rep(c(attr(terms, "term.labels")), tableTerms)
  rowIndex <- list()
  for(i in 1:length(c(tt$between, tt$within))){
    rowIndex[[i]] <- which(termsLabels == c(tt$between, tt$within)[i])
    names(rowIndex)[i] <- c(tt$between, tt$within)[i]
  }
  return(lapply(rowIndex, function(x){
    L[x,, drop = FALSE]
    }))
}

#' Refit Model with Restricted Maximum Likelihood
#'
#' @export
#'
#' @description
#' Models fit with maximum likelihood are refit with restricted maximum
#' likelihood (REML).  In this package, models are fit with REML to ensure the
#' (co)variance parameter estimates are unbiased estimates.
#'
#' @param fMod A model of class \code{\link[nlme]{gls}} or \code{\link[nlme]{lme}}.
#'
#' @returns
#' This function will return a model of class \code{\link[nlme]{gls}} or
#' \code{\link[nlme]{lme}} fit with restricted maximum likelihood.
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong, method = "ML")
#'
#' # Marginal variance-covariance matrix of the random/repeated effects for
#' # SSCP E matrix
#' Vmatrix(glsMod)$V
#'
#' # Fit with REML for marginal variance-covariance matrix of the
#' # random/repeated effects for compound symmetry covariance
#' Vmatrix(MLtoREML(update(glsMod,
#'                         correlation = corCompSymm(form = ~ 1 | Subject),
#'                         weights = NULL)))

MLtoREML <- function(fMod){
  if((!any(class(fMod) %in% c("gls", "lme"))))
    stop('The model must be a "gls" or "lme" object.', call. = FALSE)
  # Refit models w/ REML if fit by ML
  if(class(fMod) %in% c("gls", "lme")){
    if(fMod$method == "ML"){
      fMod <- update(fMod, method = "REML")
      warning("Refitting model with restricted maximum likelihood.", call. = FALSE)
    }
  }
  return(fMod)
}

#' Refit Model with Maximum Likelihood
#'
#' @export
#'
#' @description
#' Models fit with restricted maximum likelihood are refit with maximum
#' likelihood (ML).  In this package, models are fit with ML for likelihood
#' ratio tests and to calculate the error sum of squares and crossproducts
#' matrix (**E**).
#'
#' @param fMod A model of class \code{\link[nlme]{gls}} or \code{\link[nlme]{lme}}.
#'
#' @returns
#' This function will return a model of class \code{\link[nlme]{gls}} or
#' \code{\link[nlme]{lme}} fit with maximum likelihood.
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#'
#' # Likelihood ratio test
#' anova(REMLtoML(update(glsMod, model = Distance ~ Sex + Age)),
#'       REMLtoML(glsMod))
#' wilksLambda(glsMod, approximation = "LR")$withinTests$wilks[2, ]
#'
#' # SSCP E matrix
#' 27 * Vmatrix(REMLtoML(glsMod))$V
#' Ematrix(glsMod)$EB

REMLtoML <- function(fMod){
  if((!any(class(fMod) %in% c("gls", "lme"))))
    stop('The model must be a "gls" or "lme" object.', call. = FALSE)
  # Refit models w/ ML if fit by REML
  if(class(fMod) %in% c("gls", "lme")){
    if(fMod$method == "REML"){
      fMod <- update(fMod, method = "ML")
      warning("Refitting model with maximum likelihood.", call. = FALSE)
    }
  }
  return(fMod)
}

Vmatrix.ulm <- function(fMod, individual = NULL, full = FALSE){
  if(class(fMod) %in% "gls"){
    if(full){
      Vmatrices <- lapply(names(table(nlme::getGroups(fMod))), function(individuals){
        nlme::getVarCov(fMod, type = "marginal", individual = individuals)
      })
      Vmatrices <- lapply(Vmatrices, unclass)
      Vmatrix <- as.matrix(Matrix::bdiag(Vmatrices))
    }else{
      if(is.null(individual))
        individual <- names(which.max(table(nlme::getGroups(fMod))))
      Vmatrix <- nlme::getVarCov(fMod, type = "marginal", individual = individual)
      Vmatrix <- unclass(Vmatrix)
    }
  }else if(class(fMod) %in% "lme"){
    if(full){
      Vmatrices <- lapply(names(table(nlme::getGroups(fMod))), function(individuals){
        nlme::getVarCov(fMod, type = "marginal", individuals = individuals)[[1]]
      })
      Vmatrix <- as.matrix(Matrix::bdiag(Vmatrices))
    }else{
      if(is.null(individual))
        individual <- names(which.max(table(nlme::getGroups(fMod))))
      Vmatrix <- nlme::getVarCov(fMod, type = "marginal", individuals = individual)[[1]]
    }
  }

  if(mmemmuris::covStruct(fMod) == "cs" & full == FALSE){
    CS <- unique(Vmatrix[row(Vmatrix) != col(Vmatrix)])
    MSresidual <- unique(diag(Vmatrix)) - CS
    V <- list(V = Vmatrix, csCovParms = c(CS, MSresidual))
  }else{
    V <- list(V = Vmatrix)
  }
  class(V) <- c("list", "Vmatrix.ulm")
  return(V)
}

Ematrix.ulm <- function(fMod, individual = NULL, ss = NULL){
  if(mmemmuris::covStruct(fMod) %in% c("cs", "other"))
    stop("Only unstructured covariance models are allowed.", call. = FALSE)
  fMod <- mmemmuris::REMLtoML(fMod)
  Vmatrix <- mmemmuris:::Vmatrix.ulm(fMod, individual = individual)
  tt <- mmemmuris::termsType(fMod)
  data <- mmemmuris:::getDataset(fMod)
  if(!is.null(ss) & is.numeric(ss) & length(ss) == 1L)
    n <- ss
  else
    n <- mmemmuris::completeData(fMod)$n
  if(ncol(Vmatrix$V) == 0L)
    stop('Please provide a valid V matrix through the "individual" argument.', call. = FALSE)
  MB <- diag(1, nrow(Vmatrix$V))
  EB <- n * t(MB) %*% Vmatrix$V %*% MB
  M <- contr.sum(nrow(Vmatrix$V))
  E <- n * t(M) %*% Vmatrix$V %*% M
  SSE <- sum(diag(E %*% solve(t(M) %*% M)))
  Ematrix <- list(E = E, M = M, SSE = SSE, EB = EB, MB = MB, Vmatrix = Vmatrix$V,
                  X = tt$Xbetween, n = n, r = nrow(MB))
  class(Ematrix) <- c("list", "Ematrix.ulm")
  return(Ematrix)
}

Ematrix.mlm <- function(fMod){
  if(is.null(fMod$call$data))
    stop('Please use the "data" argument in the multivariate model.', call. = FALSE)
  dat <- na.omit(eval(fMod$call$data))
  cols <- colnames(coef(fMod))
  Y <- as.matrix(dat[, cols, drop = FALSE])

  formulae <- paste(attr(terms(fMod), "term.labels"), collapse = "+")
  if(formulae == ""){
    formulae <- "1"
    X <- model.matrix(as.formula(paste("~", formulae)), dat)
    B <- t(as.matrix(apply(Y, 2, function(y){
      MASS::ginv(t(X) %*% X) %*% t(X) %*% as.matrix(y)
    })))
  }else{
    X <- model.matrix(as.formula(paste("~", formulae)), dat)
    B <- apply(Y, 2, function(y){
      MASS::ginv(t(X) %*% X) %*% t(X) %*% y
    })
  }
  MB <- diag(1, length(cols))
  M <- contr.sum(length(cols))
  EB <- t(MB) %*% (t(Y) %*% Y - t(B) %*% (t(X) %*% X) %*% B) %*% MB
  E <- t(M) %*% (t(Y) %*% Y - t(B) %*% (t(X) %*% X) %*% B) %*% M
  SSE <- sum(diag(E %*% solve(t(M) %*% M)))
  Ematrix <- list(E = E, M = M, SSE = SSE, EB = EB, MB = MB, B = B, Y = Y, X = X, n = nrow(dat), r = nrow(EB))
  class(Ematrix) <- c("list", "Ematrix.mlm")
  return(Ematrix)
}

#' MANOVA/RM ANOVA SSCP Error (**E**) Matrix
#'
#' @export
#'
#' @description
#' This function returns the error (**E**) sums of squares and crossproducts (SSCP)
#' matrices from a marginal or linear mixed-effects model.  The **E** matrix is
#' calculated from the marginal variance-covariance (**V**) matrix of an
#' unstructured covariance model fit with maximum likelihood.  The **E** matrix is
#' used to calculate multivariate tests for MANOVA, univariate tests for
#' RM ANOVA, sphericity tests for within-subject effects, and epsilon
#' corrections for violations of sphericity.
#'
#' @param fMod A model of class \code{\link[nlme]{gls}}, \code{\link[nlme]{lme}}, or \code{mlm}.
#' @param individual The subject whose marginal variance-covariance matrix
#' should be used, where **V = ZGZ' + R** (**Z** is the random effects design
#' matrix, **G** is the corresponding (co)variances of the random effects matrix,
#' and **R** represents the (co)variances of the repeated effects matrix).  The
#' default is the first subject in the factor.
#'
#' @returns
#' \tabular{ll}{
#'    `E` \tab A sum of squares cross products (SSCP) error matrix for
#'    within-subject effects. \cr
#'    \tab \cr
#'    `M` \tab A sum-to-zero contrast matrix for within-subject effects. \cr
#'    \tab \cr
#'    `SSE` \tab The error sum of squares for within-subject effects calculated
#'    as \eqn{SS_E = trace(}**E**\eqn{(}**M'M**\eqn{) ^ {-1})}. \cr
#'    \tab \cr
#'    `EB` \tab A sum of squares cross products (SSCP) error matrix for
#'    between-subject effects. \cr
#'    \tab \cr
#'    `MB` \tab A sum-to-zero contrast matrix for between-subject effects. \cr
#'    \tab \cr
#'    `Vmatrix` \tab The **V** matrix from the marginal or mixed-effects model. \cr
#'    \tab \cr
#'    `X` \tab A between-subject effects design matrix. \cr
#'    \tab \cr
#'    `n` \tab The sample size. \cr
#'    \tab \cr
#'    `r`  \tab The total number of repeated measurements for data in the long format. \cr
#'    \tab \cr
#' }
#'
#' @details
#' The sum of squares cross products (SSCP) error matrices are
#' calculated as
#'
#' **E = M'** \eqn{(}**Y'Y**\eqn{-} **\eqn{\hat{B}}'** \eqn{(}**X'X**\eqn{)}**\eqn{\hat{B}}** \eqn{)}**M**
#'
#' where **M** is an identity matrix for between-subject effects and a
#' sum-to-zero contrast matrix for within-subject effects, **Y** is a combined
#' matrix of observations, **\eqn{\hat{B}}** is a matrix of combined coefficients from the
#' ordinary least squares models, and **X** is the design matrix of the
#' independent variables.
#'
#' For marginal and mixed-effects models, the SSCP **E** matrix is defined as
#'
#' **E =** \eqn{n}**M'VM**
#'
#' where \eqn{n} is the number of subjects, **M** is an identity matrix for
#' between-subject effects and a sum-to-zero contrast matrix for within-subject
#' effects, and **V** is the
#' marginal variance-covariance matrix where
#'
#' **V = ZGZ' + R** (**Z** is the random effects design matrix, **G** is the
#' corresponding (co)variances of the random effects matrix, and **R**
#' represents the (co)variances of the repeated effects matrix) for a
#' \code{\link[nlme]{gls}} or \code{\link[nlme]{lme}}
#' unstructured covariance model object (refit with maximum likelihood instead
#' of restricted maximum likelihood).
#'
#' @references {\url{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.4/statug/statug_introreg_sect038.htm}}
#' @references {\url{https://www.lesahoffman.com/PSYC943/mv12psyc943_lecture13.pdf}}
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#'
#' # Multivariate linear model
#' manovaMod <- manova(cbind(Distance8, Distance10, Distance12, Distance14) ~ Sex,
#'                     data = orthodontWide)
#'
#' # SSCP between-subject E matrix
#' Ematrix(glsMod)$EB
#'
#' # SSCP between-subject E matrix
#' Ematrix(manovaMod)$EB
#'
#' # SSCP within-subject E matrix
#' Ematrix(glsMod)$E
#'
#' # SSCP within-subject E matrix
#' Ematrix(manovaMod)$E
#'
#' @seealso
#' {\code{\link[mmemmuris]{coefs}}, \code{\link[mmemmuris]{Vmatrix}}}

Ematrix <- function(fMod, individual = NULL, ss = NULL){
  if(any(class(fMod) %in% c("gls", "lme"))){
    Ematrices <- mmemmuris:::Ematrix.ulm(fMod, individual, ss)
  }else if(any(class(fMod) %in% "mlm")){
    Ematrices <- mmemmuris:::Ematrix.mlm(fMod)
  }else
    stop('Please provide a model of class "gls" or "lme".', call. = FALSE)
  return(Ematrices)
}

Hmatrix.mlm <- function(fMod, within = "Time"){
  if(!(any(class(fMod) %in% "mlm")))
    stop("Please provide a multivariate linear model.", call. = FALSE)
  if(is.null(fMod$call$data))
    stop('Please use the "data" argument in the multivariate model.', call. = FALSE)
  dat <- na.omit(eval(fMod$call$data))
  cols <- colnames(coef(fMod))
  Y <- as.matrix(dat[, cols])

  formulae <- paste(attr(terms(fMod), "term.labels"), collapse = "+")
  if(formulae == ""){
    formulae <- "1"
    X <- model.matrix(as.formula(paste("~", formulae)), dat)
    B <- t(as.matrix(apply(Y, 2, function(y){
      MASS::ginv(t(X) %*% X) %*% t(X) %*% as.matrix(y)
    })))
  }else{
    X <- model.matrix(as.formula(paste("~", formulae)), dat)
    B <- apply(Y, 2, function(y){
      MASS::ginv(t(X) %*% X) %*% t(X) %*% y
    })
  }

  Lbetween <- mmemmuris:::Lmatrix.mlm(fMod)

  nr <- nrow(t(B))

  MB <- diag(1, nr)
  M <- contr.sum(nr)

  HB <- lapply(Lbetween, function(x){
    t(MB) %*% (t(B) %*% t(x)) %*% solve(x %*% MASS::ginv(t(X) %*% X) %*% t(x)) %*% x %*% (B) %*% MB
  })

  H <- lapply(Lbetween, function(x){
    t(M) %*% (t(B) %*% t(x)) %*% solve(x %*% MASS::ginv(t(X) %*% X) %*% t(x)) %*% x %*% (B) %*% M
  })
  names(H) <- paste(names(H), within, sep = ":")
  names(H) <- gsub("\\(Intercept\\):", "", names(H))

  SSH <- sapply(H, function(x){
    sum(diag(x %*% solve(t(M) %*% M)))
  })

  Hmatrix <- list(H = H, M = M, SSH = SSH, HB = HB, MB = MB, B = B, Y = Y, X = X,
                  Lbetween = Lbetween, n = nrow(dat))
  class(Hmatrix) <- c("list", "Hmatrix.mlm")
  return(Hmatrix)
}

Hmatrix.ulm <- function(fMod){
  tt <- mmemmuris::termsType(fMod)
  dataset <- na.omit(mmemmuris:::getDataset(fMod))
  Xbetween <- tt$Xbetween
  X <- model.matrix(fMod, data = dataset)
  B <- list()
  orders <- attr(terms(fMod), "order")
  names(orders) <- attr(terms(fMod), "term.labels")
  withinFac <- names(orders[(names(orders) %in% tt$within) & orders ==
                              1])
  withinLevels <- droplevels(dataset[[withinFac]])
  for(i in 1:length(table(withinLevels))){
    dataset[[withinFac]] <- relevel(withinLevels, ref = i)
    fMod <- update(fMod, data = dataset)
    B[[i]] <- mmemmuris::coefs(fMod)
  }
  B <- do.call("cbind", B)

  Lbetween <- t(as.matrix(Matrix::expand(Matrix::lu(t(Xbetween) %*% Xbetween))$L))
  rownames(Lbetween) <- colnames(Lbetween) <- colnames(Xbetween)

  termsIndex <- attr(Xbetween, "assign")
  tableTerms <- table(termsIndex)
  termsLabels <- rep(c("(Intercept)", attr(terms(fMod), "term.labels")[attr(terms(fMod), "term.labels") %in% tt$between]), tableTerms)
  term <- c("(Intercept)", attr(terms(fMod), "term.labels")[attr(terms(fMod), "term.labels") %in% tt$between])
  rowIndex <- list()
  for(i in 1:length(term)){
    rowIndex[[i]] <- which(termsLabels == term[i])
    names(rowIndex)[i] <- term[i]
  }
  Lbetween <- lapply(rowIndex, function(x){
    Lbetween[x,, drop = FALSE]
  })

  termsIndex2 <- attr(X, "assign")

  tableTerms2 <- table(termsIndex2)

  termsLabels2 <- rep(c("(Intercept)", attr(terms(fMod), "term.labels")), tableTerms2)
  term <- c("(Intercept)", attr(terms(fMod), "term.labels")[attr(terms(fMod), "term.labels") %in% tt$between])

  rowIndex <- list()
  for(i in 1:length(term)){
    rowIndex[[i]] <- which(termsLabels2 == term[i])
    names(rowIndex)[i] <- term[i]
  }

  B <- B[unlist(rowIndex),, drop = FALSE]

  nr <- nrow(t(B))

  MB <- diag(1, nr)
  M <- contr.sum(nr)

  H <- lapply(Lbetween, function(x){
    t(M) %*% t(B) %*% t(x) %*% solve(x %*% MASS::ginv(t(Xbetween) %*% Xbetween) %*% t(x)) %*% x %*% B %*% M
  })

  H <- H[1:length(tt$within)]

  SSH <- sapply(H, function(x){
    sum(diag(x %*% solve(t(M) %*% M)))
  })

  HB <- lapply(Lbetween, function(x){
    t(MB) %*% t(B) %*% t(x) %*% solve(x %*% MASS::ginv(t(Xbetween) %*% Xbetween) %*% t(x)) %*% x %*% B %*% MB
  })

  names(H) <- tt$within

  Hmatrix <- list(H = H, M = M, SSH = SSH, HB = HB, MB = MB, B = B, X = Xbetween,
                  Lbetween = Lbetween)
  class(Hmatrix) <- c("list", "Hmatrix.ulm")
  return(Hmatrix)
}

#' MANOVA/RM ANOVA SSCP Hypothesis (**H**) Matrix
#' @export
#'
#' @description
#' This function returns the hypothesis (**H**) sums of squares and crossproducts (SSCP)
#' matrices from a marginal or linear mixed-effects model.  The **H** matrix is
#' used to calculate multivariate tests for MANOVA and univariate tests for
#' RM ANOVA.
#'
#' @param fMod A model of class \code{\link[nlme]{gls}}, \code{\link[nlme]{lme}}, or \code{mlm}.
#' @param within The name that should be used for the within-subject effect.
#' The default is "Time".
#'
#' @details
#' The sum of squares and crossproducts (SSCP) hypothesis matrix is
#'
#' **H** \eqn{=} **M'** \eqn{(}**L\eqn{\hat{B}}** \eqn{)}'\eqn{(}**L**\eqn{(}**X'X**\eqn{) ^ {-}}**L'** \eqn{) ^ {-1}}**L\eqn{\hat{B}}M**
#'
#' where **M** is a sum-to-zero contrast matrix, **Y** is a combined
#' matrix of observations, **\eqn{\hat{B}}** is a matrix of combined coefficients from the
#' ordinary least squares models, **X** is the design matrix of the
#' independent variables, and the contrast coefficients matrix (**L**) is
#' calculated from the **LU** decomposition of the crossproducts matrix
#' (**X'X**).
#'
#' For marginal and mixed effects models, the **H** matrix is calculated with
#' the same formula as above.  The **\eqn{\hat{B}}** matrix is calculated by switching the
#' reference levels of the within-subject effect, refitting the models, and
#' extracting the model coefficients.
#'
#' @returns
#' \tabular{ll}{
#'    `H`  \tab The sums of squares and crossproducts hypothesis matrices for within-subject effects. \cr
#'    `M` \tab A sum-to-zero contrast matrix for within-subject effects. \cr
#'    `SSH` \tab The sums of squares of the hypothesis matrices for within-subject effects. \cr
#'    `HB`  \tab The sums of squares and crossproducts hypothesis matrices for between-subject effects. \cr
#'    `MB` \tab An identity contrast matrix for between-subject effects. \cr
#'    `B` \tab The multivariate linear model coefficients. \cr
#'    `X` \tab The design matrix of the between-subject effects. \cr
#'    `L`  \tab The contrast coefficients matrices. \cr
#' \tab \cr
#' }
#'
#' @references {\url{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.4/statug/statug_introreg_sect038.htm}}
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#'
#' # Multivariate linear model
#' manovaMod <- manova(cbind(Distance8, Distance10, Distance12, Distance14) ~ Sex,
#'                     data = orthodontWide)
#'
#' # SSCP between-subject H matrix
#' Hmatrix(glsMod)$HB
#'
#' # SSCP between-subject H matrix
#' Hmatrix(manovaMod)$HB
#'
#' # SSCP within-subject H matrix
#' Hmatrix(glsMod)$H
#'
#' # SSCP within-subject H matrix
#' Hmatrix(manovaMod)$H

Hmatrix <- function(fMod, within = "Time"){
  if(any(class(fMod) %in% c("gls", "lme"))){
    Hmatrices <- mmemmuris:::Hmatrix.ulm(fMod)
  }else if(any(class(fMod) %in% "mlm")){
    Hmatrices <- mmemmuris:::Hmatrix.mlm(fMod, within)
  }else
    stop('Please provide a model of class "gls", "lme", or "mlm".', call. = FALSE)
  return(Hmatrices)
}

#' Get Coefficients from Model
#'
#' @export
#'
#' @description
#' This function will return the (generalized) least squares estimates of
#' parameters in marginal and mixed-effects models and multivariate linear
#' models.
#'
#' @details
#' Marginal and mixed-effects model (fixed-effects) coefficients can be
#' calculated with generalized least squares as
#'
#' **\eqn{\hat{\beta}}** \eqn{ = (}**X'V** \eqn{^ {-1}}**X**\eqn{) ^ {-}}**X'V** \eqn{^ {-1}}**y**
#'
#' where **X** is the design matrix of the independent variables,
#' **V** is the marginal variance-covariance matrix of the random/repeated
#' effects where **V = ZGZ' + R** (**Z** is the random effects design matrix,
#' **G** is the corresponding (co)variances of the random effects matrix, and
#' **R** represents the (co)variances of the repeated effects matrix), and **y**
#' is a vector of observations.
#'
#' Multivariate linear model coefficients can be calculated with ordinary least
#' squares as
#'
#' **\eqn{\hat{B}}** \eqn{ = (}**X_b'X_b** \eqn{) ^ {-}}**X_b'Y**
#'
#' where **X_b** is the design matrix of the between-subject effects and
#' **Y** is a combined matrix of observations.
#'
#' To get the multivariate linear model coefficients from a \code{\link[nlme]{gls}} or \code{\link[nlme]{lme}}
#' model, use the \link[mmemmuris:Hmatrix]{Hmatrix} function as it calls this
#' function for each level of the within-subject effect.
#'
#' @param fMod A model of class \code{\link[nlme]{gls}}, \code{\link[nlme]{lme}},
#' or \code{mlm}.
#'
#' @returns
#' This function will return a matrix of parameter estimates for the model.
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#' # Model coefficients
#' coefs(glsMod)
#'
#' # Linear mixed-effects model
#' # Unstructured covariance
#' lmeMod <- lme(Distance ~ Sex * Age,
#'               random = ~ Age | Subject,
#'               na.action = na.omit,
#'               data = orthodontLong)
#' # Fixed-effects model coefficients
#' coefs(lmeMod)
#'
#' # Multivariate linear model
#' manovaMod <- manova(cbind(Distance8, Distance10, Distance12, Distance14) ~ Sex,
#'                     data = orthodontWide)
#' # Model coefficients
#' coefs(manovaMod)
#'
#' # Get multivariate linear model coefficients from a univariate model
#' coefs(glsMod)[1:2, ]
#' orthodontLong$Age <- relevel(orthodontLong$Age, ref = "10")
#' coefs(update(glsMod, data = orthodontLong))[1:2, ]
#' orthodontLong$Age <- relevel(orthodontLong$Age, ref = "12")
#' coefs(update(glsMod, data = orthodontLong))[1:2, ]
#' orthodontLong$Age <- relevel(orthodontLong$Age, ref = "14")
#' coefs(update(glsMod, data = orthodontLong))[1:2, ]
#'
#' @seealso
#' {\code{\link[mmemmuris]{Hmatrix}}, \code{\link[mmemmuris]{Ematrix}}, \code{\link[mmemmuris]{waldF}}, \code{\link[mmemmuris]{uniTest}}}

coefs <- function(fMod){
  if(any(class(fMod) %in% c("gls", "lme"))){
    if(class(fMod) %in% "lme")
      coefs <- cbind(nlme::fixef(fMod))
    else
      coefs <- cbind(coef(fMod))
  }else if(any(class(fMod) %in% "mlm"))
    coefs <- coef(fMod)
  else
    stop('The model must be a "gls", "lme", or "mlm" object.', call. = FALSE)
  return(coefs)
}

#' Mauchly's Sphericity Tests for Within-Subject Effects
#'
#' @export
#'
#' @description
#' This function will calculate Mauchly's test of sphericity for within-subject
#' effects for marginal models (\code{\link[nlme]{gls}}), mixed-effects models
#' (\code{\link[nlme]{lme}}), and multivariate linear models (\code{mlm}).
#'
#' @param Ematrices An object of class `Ematrix.ulm` or `Ematrix.mlm`.
#'
#' @details
#' The sum of squares cross products (SSCP) error matrices are
#' calculated as
#'
#' **E = M'** \eqn{(}**Y'Y**\eqn{-} **\eqn{\hat{B}}'** \eqn{(}**X'X**\eqn{)}**\eqn{\hat{B}}** \eqn{)}**M**
#'
#' where **M** is an identity matrix for between-subject effects and a
#' sum-to-zero contrast matrix for within-subject effects, **Y** is a combined
#' matrix of observations, **\eqn{\hat{B}}** is a matrix of combined coefficients from the
#' ordinary least squares models, and **X** is the design matrix of the
#' independent variables.
#'
#' For marginal and mixed-effects models, the SSCP **E** matrix is defined as
#'
#' **E =** \eqn{n}**M'VM**
#'
#' where \eqn{n} is the number of subjects, **M** is an identity matrix for
#' between-subject effects and a sum-to-zero contrast matrix for within-subject
#' effects, and **V** is the
#' marginal variance-covariance matrix where
#'
#' **V = ZGZ' + R** (**Z** is the random effects design matrix, **G** is the
#' corresponding (co)variances of the random effects matrix, and **R**
#' represents the (co)variances of the repeated effects matrix) for a
#' \code{\link[nlme]{gls}} or \code{\link[nlme]{lme}}
#' unstructured covariance model object (refit with maximum likelihood instead
#' of restricted maximum likelihood).
#'
#' The within-subject SSCP **E** matrix is used to calculate Mauchly's test for
#' sphericity for within-subject effects.
#'
#' Mauchly's test for sphericity is calculated as
#'
#' \eqn{W = (\Pi\lambda_i of} **E**\eqn{(}**M'M**\eqn{) ^ {-1}) / (((\Sigma\lambda_i of} **E**\eqn{(}**M'M**\eqn{) ^ {-1}) / (r - 1)) ^ (r - 1))}
#'
#' Let \eqn{r} be equal to the number of repeated measurements, \eqn{v}
#' be the residual degrees of freedom (between-subject design matrix), and the correction factor equal
#' \deqn{\rho = 1 - (2(r - 1) ^ 2 + r + 1) / (6(r - 1)v)}  The \eqn{\chi ^ 2} statistic is
#' equal to the product \eqn{-\rho v log(W)} with degrees of freedom
#' \eqn{(r(r - 1) / 2) - 1}.
#'
#' @returns
#' \tabular{ll}{
#'    `sphericityTests`  \tab A `data.frame` containing the Mauchly's sphericity
#'    test statistics. \cr
#'    `r`  \tab The total number of repeated measurements for data in the long format. \cr
#'    `v`  \tab The residual degrees of freedom. \cr
#'    `E`  \tab The sums of squares and crossproducts error matrix for within-subject effects. \cr
#'    `M`  \tab A sum-to-zero contrast matrix for within-subject effects. \cr
#' \tab \cr
#' }
#'
#' @references {\url{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.4/statug/statug_introreg_sect038.htm}}
#' @references {\url{https://www.lesahoffman.com/PSYC943/mv12psyc943_lecture13.pdf}}
#' @references {\url{https://support.sas.com/rnd/app/stat/papers/mixedglm.pdf}}
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#' sphericityTests(Ematrix(glsMod))
#'
#' # Linear mixed-effects model
#' # Unstructured covariance
#' lmeMod <- lme(Distance ~ Sex * Age,
#'               random = ~ Age | Subject,
#'               na.action = na.omit,
#'               data = orthodontLong)
#' sphericityTests(Ematrix(lmeMod))
#'
#' # Multivariate linear model
#' manovaMod <-
#'   manova(cbind(Distance8, Distance10, Distance12, Distance14) ~ Sex,
#'          data = orthodontWide)
#' sphericityTests(Ematrix(manovaMod))
#'
#' @seealso
#' {\code{\link[mmemmuris]{coefs}}, \code{\link[mmemmuris]{Vmatrix}}, \code{\link[mmemmuris]{Ematrix}}, \code{\link[mmemmuris]{epsilon}}}

sphericityTests <- function(Ematrices){
  if(any(class(Ematrices) %in% c("Ematrix.mlm", "Ematrix.ulm"))){
    if(is.null(mmemmuris:::inverseMatrix(Ematrices$E))){
      warning("The SSCP E matrix is singular.  Sphericity tests are not available.", call. = FALSE)
      sphericityTests <- list(r = Ematrices$r)
      class(sphericityTests) <- c("list", "sphericityTests")
      return(sphericityTests)
    }
    if(1 / (Ematrices$r - 1) == 1L){
      warning("The lower-bound epsilon is 1.  Sphericity is automatically satisfied.", call. = FALSE)
      sphericityTests <- list(r = Ematrices$r)
      class(sphericityTests) <- c("list", "sphericityTests")
      return(sphericityTests)
    }
    W <- det(diag(eigen(Ematrices$E %*% solve(t(Ematrices$M) %*% Ematrices$M))$values)) / (sum(eigen(Ematrices$E %*% solve(t(Ematrices$M) %*% Ematrices$M))$values) / (Ematrices$r - 1)) ^ (Ematrices$r - 1)
    v <- Ematrices$n - c(Matrix::rankMatrix(Ematrices$X))
    rho <- 1 - (2 * (Ematrices$r - 1) ^ 2 + (Ematrices$r - 1) + 2) / (6 * (Ematrices$r - 1) * v)
    ChiW <- -rho * v * log(W)
    # Construct the results data.frame
    spher <- data.frame(W = W,
                        rho = rho,
                        df = ((Ematrices$r - 1) * (Ematrices$r) / 2) - 1,
                        Chisq = ChiW,
                        p.value = 1 - pchisq(ChiW, ((Ematrices$r - 1) * (Ematrices$r) / 2) - 1))
    rownames(spher) <- ""
    colnames(spher) <- c("W", "rho", "df", "Chisq", "Pr(>Chisq)")
    spher <- cbind(spher, mmemmuris::sigStars(spher$`Pr(>Chisq)`))
    colnames(spher)[length(colnames(spher))] <- rownames(spher)[length(rownames(spher))] <- ""
    sphericityTests <- list(sphericityTests = spher, r = Ematrices$r, v = v, E = Ematrices$E, M = Ematrices$M)
    class(sphericityTests) <- c("list", "sphericityTests")
    return(sphericityTests)
  }else
    stop('Please provide an object of class "Ematrix.mlm" or "Ematrix.ulm".', call. = FALSE)
}

#' Print a Sphericity Tests Object
#'
#' @exportS3Method
print.sphericityTests <- function(sphericityTests, digits = 4){
  if(length(sphericityTests) > 0L){
    if((sphericityTests$r - 1) > 1L){
      if(digits < 3)
        stop("No less than 3 digits can be printed.", call. = FALSE)
      cat("\nMauchly Sphericity Tests")
      cat("\n------------------------\n")
      spher <- cbind(format(round(sphericityTests$sphericityTests[, -ncol(sphericityTests$sphericityTests)], digits), nsmall = digits,
                                                      scientific = FALSE),
                                               sphericityTests$sphericityTests[, ncol(sphericityTests$sphericityTests)])
      colnames(spher)[length(colnames(spher))] <- rownames(spher) <- ""
      print(spher)
      cat("---")
      cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
    }
  }
}

colSpaceX2inColSpaceX1 <- function(X1, L){
  if(is.matrix(L)){
    if(ncol(X1) != ncol(L))
      stop("The number of columns of the design matrix and L matrix don't match.")
    else{
      tL <- t(L)
      orthComp <- qr.Q(qr(cbind(tL)), complete = TRUE)[, -c(1:Matrix::rankMatrix(tL)), drop = FALSE]
      X2 <- X1 %*% orthComp
      colnames(X2) <- paste0("Col", 1:ncol(X2))
    }
    return(X2)
  }else
    stop('"L" must be a matrix.', call. = FALSE)
}

Lmatrix.mlm <- function(fMod){
  terms <- terms(fMod)
  X <- model.matrix(fMod)
  if (ncol(X) == 1L)
    L <- matrix(1L)
  else
    L <- t(as.matrix(Matrix::expand(Matrix::lu(t(X) %*% X))$L))
  rownames(L) <- colnames(L) <- colnames(X)
  termsIndex <- attr(X, "assign")
  tableTerms <- table(termsIndex)
  if(names(tableTerms)[1] == "0"){
    termsLabels <- rep(c("(Intercept)", attr(terms, "term.labels")), tableTerms)
    term <- c("(Intercept)", attr(terms, "term.labels"))
  }else{
    termsLabels <- rep(c(attr(terms, "term.labels")), tableTerms)
    term <- attr(terms, "term.labels")
  }
  rowIndex <- list()
  for(i in 1:length(term)){
    rowIndex[[i]] <- which(termsLabels == term[i])
    names(rowIndex)[i] <- term[i]
  }
  return(lapply(rowIndex, function(x){
    L[x,, drop = FALSE]
  }))
}

Vmatrix.mlm <- function(fMod, within = "Time"){
  Ematrices <- mmemmuris:::Ematrix.mlm(fMod)

  unV <- Ematrices$EB / Ematrices$n

  dat <- na.omit(eval(fMod$call$data))
  formulae <- paste(attr(terms(fMod), "term.labels"), collapse = "+")
  if(formulae == "")
    formulae <- "1"
  dat$y <- rowSums(Ematrices$Y) / sqrt(ncol(Ematrices$Y))
  test <- aov(as.formula(paste("y ~", formulae)), data = dat)
  betweenTests <- as.data.frame(summary(test)[[1]])

  SSE <- sum(diag(Ematrices$E %*% solve(t(Ematrices$M) %*% Ematrices$M)))

  MSbetween <- betweenTests$`Mean Sq`[length(betweenTests$`Mean Sq`)]
  MSresidual <- SSE / (fMod$df.residual * ncol(Ematrices$M))

  CS <- (MSbetween - MSresidual) / nrow(Ematrices$M)

  csV <- matrix(CS, nrow(Ematrices$M), nrow(Ematrices$M)) + diag(MSresidual, nrow(Ematrices$M))

  Vmatrix <- list(unV = unV, EB = Ematrices$EB, csV = csV, csCovParms = c(CS, MSresidual), E = Ematrices$E)
  class(Vmatrix) <- c("list", "Vmatrix.mlm")
  return(Vmatrix)
}

#' Marginal Variance-Covariance Matrix for Random/Repeated Effects
#'
#' @export
#'
#' @description
#' This function calculates the marginal variance-covariance matrix of the
#' random/repeated effects for a marginal model (`gls`), mixed-effects model
#' (`lme`), and multivariate linear model (`mlm`).
#'
#' @param fMod A model of class \code{\link[nlme]{gls}}, \code{\link[nlme]{lme}}, or \code{mlm}.
#' @param individual The subject whose marginal variance-covariance matrix
#' should be used, where **V = ZGZ' + R** (**Z** is the random effects design
#' matrix, **G** is the corresponding (co)variances of the random effects matrix,
#' and **R** represents the (co)variances of the repeated effects matrix).  The
#' default is the first subject in the factor.
#' @param full If `TRUE`, print the full block diagonal matrix with all
#' experimental units.
#' @param within The name that should be used for the within-subject effect.
#' The default is "Time".
#'
#' @details
#' The marginal variance-covariance matrix of the random/repeated effects is
#' denoted as
#'
#' **V = ZGZ' + R**
#'
#' where **Z** is the random effects design matrix, **G** is the corresponding
#' (co)variances of the random effects matrix, and **R** represents the
#' (co)variances of the repeated effects matrix.
#'
#' @references {\url{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.3/statug/statug_mixed_details01.htm#statug_mixed009841}}
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsModUN <- gls(Distance ~ Sex * Age,
#'                 correlation =
#'                   corSymm(form = ~ as.numeric(Age) | Subject),
#'                 weights = varIdent(form = ~ 1 | Age),
#'                 na.action = na.omit,
#'                 data = orthodontLong, method = "ML")
#'
#' # Marginal variance-covariance matrix of the random/repeated effects for
#' # SSCP E matrix
#' Vmatrix(glsModUN)$V
#'
#' # Fit with REML for marginal variance-covariance matrix of the
#' # random/repeated effects for compound symmetry covariance
#' glsModCS <- gls(Distance ~ Sex * Age,
#'                 correlation =
#'                   corCompSymm(form = ~ 1 | Subject),
#'                 na.action = na.omit,
#'                 data = orthodontLong)
#' Vmatrix(glsModCS)
#'
#' # Multivariate linear model
#' manovaMod <- manova(cbind(Distance8, Distance10, Distance12, Distance14) ~ Sex,
#'                     data = orthodontWide)
#' Vmatrix(manovaMod)[c("unV", "csV")]

Vmatrix <- function(fMod, individual = NULL, full = FALSE, within = "Time"){
  if(any(class(fMod) %in% c("gls", "lme"))){
    Vmatrices <- mmemmuris:::Vmatrix.ulm(fMod, individual, full)
  }else if(any(class(fMod) %in% "mlm")){
    Vmatrices <- mmemmuris:::Vmatrix.mlm(fMod, within)
  }else
    stop('Please provide a model of class "gls", "lme", or "mlm".', call. = FALSE)
  return(Vmatrices)
}

#' Covariance Structure
#'
#' @export
#'
#' @description
#' This function will determine whether a compound symmetry, unstructured, or
#' other covariance structure was used for the random/repeated effects for
#' objects of class \code{\link[nlme]{gls}} or \code{\link[nlme]{lme}}.
#'
#' @details
#' This package uses two different covariance structures: compound symmetry
#' (repeated measures ANOVA) and unstructured (multivariate analysis of
#' variance).  The compound symmetry covariance structure uses two (co)variance
#' parameters: a residual error variance (\eqn{\sigma ^ 2}) and covariance
#' parameter between repeated measurements (\eqn{\sigma_r ^ 2}).  The
#' unstructured covariance structure has the maximum number possible of
#' (co)variance parameters: \deqn{r(r + 1) / 2}
#'
#' @param fMod A model of class \code{\link[nlme]{gls}} or \code{\link[nlme]{lme}}.
#'
#' @returns
#' This function will return a character saying which covariance structure the
#' model has:  \tabular{lll}{\cr
#' \code{un} \tab - \tab unstructured \cr
#' \code{cs} \tab - \tab compound symmetry \cr
#' \code{other} \tab - \tab not unstructured or compound symmetry \cr
#' }
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#' covStruct(glsMod)
#' # Compound symmetry covariance
#' glsModCS <- gls(Distance ~ Sex * Age,
#'                 correlation = corCompSymm(form = ~ 1 | Subject),
#'                 na.action = na.omit,
#'                 data = orthodontLong)
#'covStruct(glsModCS)
#'
#' # Linear mixed-effects model
#' # Unstructured covariance
#' lmeMod <- lme(Distance ~ Sex * Age,
#'               random = ~ Age | Subject,
#'               na.action = na.omit,
#'               data = orthodontLong)
#' covStruct(lmeMod)
#' # Compound symmetry covariance (random intercepts)
#' lmeModCS <- lme(Distance ~ Sex * Age,
#'                 random = ~ 1 | Subject,
#'                 na.action = na.omit,
#'                 data = orthodontLong)
#' covStruct(lmeModCS)
#'
#' @references {\url{https://support.sas.com/resources/papers/proceedings/proceedings/sugi30/198-30.pdf}}

covStruct <- function(fMod){
  mmemmuris:::modelCheck(fMod)
  if(class(fMod) %in% "gls"){
    if(("corSymm" %in% class(fMod$modelStruct$corStruct)) & !is.null(fMod$modelStruct$varStruct) & all(all.vars(attr(fMod$modelStruct$varStruct, "formula")) %in% all.vars(attr(fMod$modelStruct$corStruct, "formula"))))
      return("un")
    if(("corCompSymm" %in% class(fMod$modelStruct$corStruct)) & is.null(fMod$modelStruct$varStruct))
      return("cs")
    return("other")
  }else if(class(fMod) %in% "lme"){
    if(("corSymm" %in% class(fMod$modelStruct$corStruct)) & !is.null(fMod$modelStruct$varStruct) & all(all.vars(attr(fMod$modelStruct$varStruct, "formula")) %in% all.vars(attr(fMod$modelStruct$corStruct, "formula"))))
      return("un")
    else if(length(all.vars(attr(fMod$modelStruct$reStruct[[attr(nlme::getGroups(fMod), "label")]], "formula"))) > 0)
      return("un")
    else if(length(all.vars(attr(fMod$modelStruct$reStruct[[attr(nlme::getGroups(fMod), "label")]], "formula"))) == 0 & ("corCompSymm" %in% class(fMod$modelStruct$corStruct) | is.null(fMod$modelStruct$corStruct)) & is.null(fMod$modelStruct$varStruct))
      return("cs")
    else
      return("other")
  }
}

raoF <- function(lambda){
  if(any(class(lambda) %in% c("lambda.mlm", "lambda.ulm"))){
    v <- lambda$v

    if(!is.null(lambda$lambdaBetween)){
      rBetween <- v - ((lambda$pBetween - lambda$qBetween + 1) / 2)
      uBetween <- (lambda$pBetween * lambda$qBetween - 2) / 4
      tBetween <- ifelse(lambda$pBetween ^ 2 + lambda$qBetween ^ 2 - 5 > 0, sqrt((lambda$pBetween ^ 2 * lambda$qBetween ^ 2 - 4) / (lambda$pBetween ^ 2 + lambda$qBetween ^ 2 - 5)), 1)
      df1Between <- lambda$pBetween * lambda$qBetween
      df2Between <- rBetween * tBetween - 2 * uBetween
      FlambdaBetween <- ((1 - lambda$lambdaBetween ^ (1 / tBetween)) / (lambda$lambdaBetween ^ (1 / tBetween))) * (df2Between / df1Between)

      termsWilksBetween <- data.frame(df1 = df1Between, df2 = df2Between,
                                      lambda = lambda$lambdaBetween,
                                      Fstat = FlambdaBetween,
                                      p.value = 1 - pf(FlambdaBetween, df1Between, df2Between))
      parmsWilksBetween <- data.frame(p = lambda$pBetween, q = lambda$qBetween,
                                      s = lambda$sBetween,
                                      v = rep(v, length(lambda$pBetween)),
                                      r = rBetween, u = uBetween, t = tBetween)

      colnames(termsWilksBetween) <- c("Num df", "Den df", "lambda", "F", "Pr(>F)")

      termsWilksBetween <- cbind(termsWilksBetween, mmemmuris::sigStars(termsWilksBetween$`Pr(>F)`))

      colnames(termsWilksBetween)[length(colnames(termsWilksBetween))] <- ""

      wilksBetween <- list(wilks = termsWilksBetween[rownames(termsWilksBetween) != "(Intercept)", ],
                           parmsWilks = parmsWilksBetween[rownames(parmsWilksBetween) != "(Intercept)", ],
                           approximation = "Rao",
                           EB = lambda$EB, HB = lambda$HB,
                           Lbetween = lambda$Lbetween, L = lambda$L, type = lambda$type,
                           Vmatrix = lambda$EB / lambda$n)
    }else{
      termsWilksBetween <- parmsWilksBetween <- wilksBetween <- NULL
    }

    rWithin <- v - ((lambda$pWithin - lambda$qWithin + 1) / 2)
    uWithin <- (lambda$pWithin * lambda$qWithin - 2) / 4
    tWithin <- ifelse(lambda$pWithin ^ 2 + lambda$qWithin ^ 2 - 5 > 0, sqrt((lambda$pWithin ^ 2 * lambda$qWithin ^ 2 - 4) / (lambda$pWithin ^ 2 + lambda$qWithin ^ 2 - 5)), 1)

    df1Within <- lambda$pWithin * lambda$qWithin
    df2Within <- rWithin * tWithin - 2 * uWithin
    FlambdaWithin <- ((1 - lambda$lambdaWithin ^ (1 / tWithin)) / (lambda$lambdaWithin ^ (1 / tWithin))) * (df2Within / df1Within)

    termsWilksWithin <- data.frame(df1 = df1Within, df2 = df2Within,
                                   lambda = lambda$lambdaWithin,
                                   Fstat = FlambdaWithin,
                                   p.value = 1 - pf(FlambdaWithin, df1Within, df2Within))
    parmsWilksWithin <- data.frame(p = lambda$pWithin, q = lambda$qWithin,
                                   s = lambda$sWithin,
                                   v = rep(v, length(lambda$pWithin)),
                                   r = rWithin, u = uWithin, t = tWithin)

    colnames(termsWilksWithin) <- c("Num df", "Den df", "lambda", "F", "Pr(>F)")

    termsWilksWithin <- cbind(termsWilksWithin, mmemmuris::sigStars(termsWilksWithin$`Pr(>F)`))

    colnames(termsWilksWithin)[length(colnames(termsWilksWithin))] <- ""

    wilksWithin <- list(wilks = termsWilksWithin,
                        parmsWilks = parmsWilksWithin,
                        approximation = "Rao",
                        E = lambda$E, H = lambda$H,
                        Lwithin = lambda$Lwithin, type = lambda$type)

    if(is.null(wilksBetween))
      wilksTests <- list(wilksWithin = wilksWithin)
    else
      wilksTests <- list(wilksBetween = wilksBetween, wilksWithin = wilksWithin)
    if(any(class(lambda) %in% "lambda.mlm"))
      class(wilksTests) <- c("list", "raoF", "mlm")
    else
      class(wilksTests) <- c("list", "raoF", "ulm")
    return(wilksTests)
  }else
    stop('Please provide an object of class "lambda.ulm" or "lambda.mlm".', call. = FALSE)
}

lrChi <- function(lambda){
  if(any(class(lambda) %in% c("lambda.mlm", "lambda.ulm"))){

    if(!is.null(lambda$lambdaBetween)){
      parmsWilksBetween <- data.frame(p = lambda$pBetween, q = lambda$qBetween,
                                      s = lambda$sBetween, v = lambda$v)
      df1Between <- lambda$pBetween * lambda$qBetween
      termsWilksBetween <- data.frame(df = df1Between, lambda = lambda$lambdaBetween,
                                      Chisq = lambda$X2Between,
                                      p.value = 1 - pchisq(lambda$X2Between, df1Between))

      colnames(termsWilksBetween) <- c("df", "lambda", "Chisq", "Pr(>Chisq)")

      termsWilksBetween <- cbind(termsWilksBetween, mmemmuris::sigStars(termsWilksBetween$`Pr(>Chisq)`))
      colnames(termsWilksBetween)[length(colnames(termsWilksBetween))] <- ""

      wilksBetween <- list(wilks = termsWilksBetween[rownames(termsWilksBetween) != "(Intercept)", ],
                           parmsWilks = parmsWilksBetween[rownames(parmsWilksBetween) != "(Intercept)", ],
                           approximation = "LR",
                           EB = lambda$EB, HB = lambda$HB,
                           Lbetween = lambda$Lbetween, L = lambda$L, type = lambda$type,
                           Vmatrix = lambda$EB / lambda$n)
    }else{
      termsWilksBetween <- parmsWilksBetween <- wilksBetween <- NULL
    }

    parmsWilksWithin <- data.frame(p = lambda$pWithin, q = lambda$qWithin,
                                   s = lambda$sWithin, v = lambda$v)
    df1Within <- lambda$pWithin * lambda$qWithin
    termsWilksWithin <- data.frame(df = df1Within, lambda = lambda$lambdaWithin,
                                   Chisq = lambda$X2Within,
                                   p.value = 1 - pchisq(lambda$X2Within, df1Within))

    colnames(termsWilksWithin) <- c("df", "lambda", "Chisq", "Pr(>Chisq)")

    termsWilksWithin <- cbind(termsWilksWithin, mmemmuris::sigStars(termsWilksWithin$`Pr(>Chisq)`))
    colnames(termsWilksWithin)[length(colnames(termsWilksWithin))] <- ""

    wilksWithin <- list(wilks = termsWilksWithin,
                        parmsWilks = parmsWilksWithin,
                        approximation = "LR",
                        E = lambda$E, H = lambda$H,
                        Lwithin = lambda$Lwithin, type = lambda$type)

    if(is.null(wilksBetween))
      wilksTests <- list(wilksWithin = wilksWithin)
    else
      wilksTests <- list(wilksBetween = wilksBetween, wilksWithin = wilksWithin)
    if(any(class(lambda) %in% "lambda.mlm"))
      class(wilksTests) <- c("list", "lrChi", "mlm")
    else
      class(wilksTests) <- c("list", "lrChi", "ulm")
    return(wilksTests)
  }else
    stop('Please provide an object of class "lambda.ulm" or "lambda.mlm".', call. = FALSE)
}

BChi <- function(lambda){
  if(any(class(lambda) %in% c("lambda.mlm", "lambda.ulm"))){

    if(!is.null(lambda$lambdaBetween)){
      parmsWilksBetween <- data.frame(p = lambda$pBetween, q = lambda$qBetween,
                                      s = lambda$sBetween, v = lambda$v)
      df1Between <- lambda$pBetween * lambda$qBetween
      X2Between <- (((lambda$pBetween - lambda$qBetween + 1) / (2)) - lambda$v) * log(lambda$lambdaBetween)
      termsWilksBetween <- data.frame(df = df1Between, lambda = lambda$lambdaBetween,
                                      Chisq = X2Between,
                                      p.value = 1 - pchisq(X2Between, df1Between))

      colnames(termsWilksBetween) <- c("df", "lambda", "Chisq", "Pr(>Chisq)")

      termsWilksBetween <- cbind(termsWilksBetween, mmemmuris::sigStars(termsWilksBetween$`Pr(>Chisq)`))
      colnames(termsWilksBetween)[length(colnames(termsWilksBetween))] <- ""

      wilksBetween <- list(wilks = termsWilksBetween[rownames(termsWilksBetween) != "(Intercept)", ],
                           parmsWilks = parmsWilksBetween[rownames(parmsWilksBetween) != "(Intercept)", ],
                           approximation = "Bartlett",
                           EB = lambda$EB, HB = lambda$HB,
                           Lbetween = lambda$Lbetween, L = lambda$L, type = lambda$type,
                           Vmatrix = lambda$EB / lambda$n)
    }else{
      termsWilksBetween <- parmsWilksBetween <- wilksBetween <- NULL
    }

    parmsWilksWithin <- data.frame(p = lambda$pWithin, q = lambda$qWithin,
                                   s = lambda$sWithin, v = lambda$v)
    df1Within <- lambda$pWithin * lambda$qWithin
    X2Within <- (((lambda$pWithin - lambda$qWithin + 1) / (2)) - lambda$v) * log(lambda$lambdaWithin)
    termsWilksWithin <- data.frame(df = df1Within, lambda = lambda$lambdaWithin,
                                   Chisq = X2Within,
                                   p.value = 1 - pchisq(X2Within, df1Within))

    colnames(termsWilksWithin) <- c("df", "lambda", "Chisq", "Pr(>Chisq)")

    termsWilksWithin <- cbind(termsWilksWithin, mmemmuris::sigStars(termsWilksWithin$`Pr(>Chisq)`))
    colnames(termsWilksWithin)[length(colnames(termsWilksWithin))] <- ""

    wilksWithin <- list(wilks = termsWilksWithin,
                        parmsWilks = parmsWilksWithin,
                        approximation = "Bartlett",
                        E = lambda$E, H = lambda$H,
                        Lwithin = lambda$Lwithin, type = lambda$type)

    if(is.null(wilksBetween))
      wilksTests <- list(wilksWithin = wilksWithin)
    else
      wilksTests <- list(wilksBetween = wilksBetween, wilksWithin = wilksWithin)
    if(any(class(lambda) %in% "lambda.mlm"))
      class(wilksTests) <- c("list", "BChi", "mlm")
    else
      class(wilksTests) <- c("list", "BChi", "ulm")
    return(wilksTests)
  }else
    stop('Please provide an object of class "lambda.ulm" or "lambda.mlm".', call. = FALSE)
}

makeWilksLambda.ulm <- function(lambda, approximation = c("Rao", "Bartlett", "LR")){
  approximation <- match.arg(approximation)

  if(approximation == "Rao"){
    raoFApprox <- raoF(lambda)

    withinTests <- raoFApprox$wilksWithin

    if(!is.null(raoFApprox$wilksBetween)){
      betweenTests <- raoFApprox$wilksBetween
      wilks <- list(betweenTests = betweenTests,
                    withinTests = withinTests,
                    LR = lambda$LR)
      class(wilks) <- c("list", "multiTest", "wilks", "ulm")
    }else{
      wilks <- list(withinTests = withinTests,
                    LR = lambda$LR)
      class(wilks) <- c("list", "multiTest", "wilks", "ulm")
    }
    return(wilks)
  }else if(approximation == "LR"){
    lrChiApprox <- lrChi(lambda)

    withinTests <- lrChiApprox$wilksWithin

    if(!is.null(lrChiApprox$wilksBetween)){
      betweenTests <- lrChiApprox$wilksBetween
      wilks <- list(betweenTests = betweenTests,
                    withinTests = withinTests,
                    LR = lambda$LR)
      class(wilks) <- c("list", "multiTest", "wilks", "ulm")
    }else{
      wilks <- list(withinTests = withinTests,
                    LR = lambda$LR)
      class(wilks) <- c("list", "multiTest", "wilks", "ulm")
    }
    return(wilks)
  }else if(approximation == "Bartlett"){
    BChiApprox <- BChi(lambda)

    withinTests <- BChiApprox$wilksWithin

    if(!is.null(BChiApprox$wilksBetween)){
      betweenTests <- BChiApprox$wilksBetween
      wilks <- list(betweenTests = betweenTests,
                    withinTests = withinTests,
                    LR = lambda$LR)
      class(wilks) <- c("list", "multiTest", "wilks", "ulm")
    }else{
      wilks <- list(withinTests = withinTests,
                    LR = lambda$LR)
      class(wilks) <- c("list", "multiTest", "wilks", "ulm")
    }
    return(wilks)
  }
}

mullerF <- function(V){
  if(any(class(V) %in% c("V.mlm", "V.ulm"))){

    if(!is.null(V$VBetween)){
      dfNumBetween <- ((V$qBetween * V$pBetween) / (V$sBetween * (V$v + V$qBetween))) * (((V$sBetween * (V$v + V$sBetween - V$pBetween) * (V$v + V$qBetween + 2) * (V$v + V$qBetween - 1)) / (V$v * (V$v + V$qBetween - V$pBetween))) - 2)
      dfDenBetween <- ((V$v + V$sBetween - V$pBetween) / (V$v + V$qBetween)) * (((V$sBetween * (V$v + V$sBetween - V$pBetween) * (V$v + V$qBetween + 2) * (V$v + V$qBetween - 1)) / (V$v * (V$v + V$qBetween - V$pBetween))) - 2)
      FVBetween <- ((dfDenBetween) / (dfNumBetween)) * ((V$VBetween) / (V$sBetween - V$VBetween))

      termsPillaiBetween <- data.frame(df1 = dfNumBetween,
                                       df2 = dfDenBetween,
                                       VBetween = V$VBetween,
                                       FVBetween = FVBetween,
                                       p.value = 1 - pf(FVBetween, dfNumBetween, dfDenBetween))
      parmsPillaiBetween <- data.frame(p = V$pBetween, q = V$qBetween,
                                       s = V$sBetween,
                                       v = V$v)

      colnames(termsPillaiBetween) <- c("Num df", "Den df", "V", "F", "Pr(>F)")

      termsPillaiBetween <- cbind(termsPillaiBetween, mmemmuris::sigStars(termsPillaiBetween$`Pr(>F)`))

      colnames(termsPillaiBetween)[length(colnames(termsPillaiBetween))] <- ""

      rownames(termsPillaiBetween) <- rownames(parmsPillaiBetween)

      if(any(class(V) %in% "V.ulm")){
        if(V$LR == TRUE){
          pillaiBetween <- list(pillai = termsPillaiBetween[(rownames(termsPillaiBetween) != "(Intercept)") & V$sBetween == 1, ],
                                parmsPillai = parmsPillaiBetween[(rownames(termsPillaiBetween) != "(Intercept)") & V$sBetween == 1, ],
                                approximation = "Muller",
                                EB = V$EB, HB = V$HB,
                                Lbetween = V$Lbetween, L = V$L, type = V$type,
                                Vmatrix = V$EB / V$n)
        }else{
          pillaiBetween <- list(pillai = termsPillaiBetween[rownames(termsPillaiBetween) != "(Intercept)", ],
                                parmsPillai = parmsPillaiBetween[rownames(parmsPillaiBetween) != "(Intercept)", ],
                                approximation = "Muller",
                                EB = V$EB, HB = V$HB,
                                Lbetween = V$Lbetween, L = V$L, type = V$type,
                                Vmatrix = V$EB / V$n)
        }
      }else{
        pillaiBetween <- list(pillai = termsPillaiBetween[rownames(termsPillaiBetween) != "(Intercept)", ],
                              parmsPillai = parmsPillaiBetween[rownames(parmsPillaiBetween) != "(Intercept)", ],
                              approximation = "Muller",
                              EB = V$EB, HB = V$HB,
                              Lbetween = V$Lbetween, L = V$L, type = V$type,
                              Vmatrix = V$EB / V$n)
      }
    }else{
      termsPillaiBetween <- parmsPillaiBetween <- pillaiBetween <- NULL
    }

    dfNumWithin <- ((V$qWithin * V$pWithin) / (V$sWithin * (V$v + V$qWithin))) * (((V$sWithin * (V$v + V$sWithin - V$pWithin) * (V$v + V$qWithin + 2) * (V$v + V$qWithin - 1)) / (V$v * (V$v + V$qWithin - V$pWithin))) - 2)
    dfDenWithin <- ((V$v + V$sWithin - V$pWithin) / (V$v + V$qWithin)) * (((V$sWithin * (V$v + V$sWithin - V$pWithin) * (V$v + V$qWithin + 2) * (V$v + V$qWithin - 1)) / (V$v * (V$v + V$qWithin - V$pWithin))) - 2)
    FVWithin <- ((dfDenWithin) / (dfNumWithin)) * ((V$VWithin) / (V$sWithin - V$VWithin))

    termsPillaiWithin <- data.frame(df1 = dfNumWithin,
                                    df2 = dfDenWithin,
                                    VWithin = V$VWithin,
                                    FVWithin = FVWithin,
                                    p.value = 1 - pf(FVWithin, dfNumWithin, dfDenWithin))
    parmsPillaiWithin <- data.frame(p = V$pWithin, q = V$qWithin,
                                    s = V$sWithin,
                                    v = V$v)

    colnames(termsPillaiWithin) <- c("Num df", "Den df", "V", "F", "Pr(>F)")

    termsPillaiWithin <- cbind(termsPillaiWithin, mmemmuris::sigStars(termsPillaiWithin$`Pr(>F)`))

    colnames(termsPillaiWithin)[length(colnames(termsPillaiWithin))] <- ""

    rownames(termsPillaiWithin) <- rownames(parmsPillaiWithin)

    if(any(class(V) %in% "V.ulm")){
      if(V$LR == TRUE){
        pillaiWithin <- list(pillai = termsPillaiWithin[V$sWithin == 1, ],
                             parmsPillai = parmsPillaiWithin,
                             approximation = "Muller",
                             E = V$E, H = V$H,
                             Lwithin = V$Lwithin, type = V$type)
      }else{
        pillaiWithin <- list(pillai = termsPillaiWithin,
                             parmsPillai = parmsPillaiWithin,
                             approximation = "Muller",
                             E = V$E, H = V$H,
                             Lwithin = V$Lwithin, type = V$type)
      }
    }else{
      pillaiWithin <- list(pillai = termsPillaiWithin,
                           parmsPillai = parmsPillaiWithin,
                           approximation = "Muller",
                           E = V$E, H = V$H,
                           Lwithin = V$Lwithin, type = V$type)
    }

    if(is.null(pillaiBetween))
      pillaiTests <- list(pillaiWithin = pillaiWithin)
    else
      pillaiTests <- list(pillaiBetween = pillaiBetween, pillaiWithin = pillaiWithin)
    if(any(class(V) %in% "V.mlm"))
      class(pillaiTests) <- c("list", "mullerF", "mlm")
    else
      class(pillaiTests) <- c("list", "mullerF", "ulm")
    return(pillaiTests)
  }else
    stop('Please provide an object of class "V.ulm" or "V.mlm".', call. = FALSE)
}

pillaiF <- function(V){
  if(any(class(V) %in% c("V.mlm", "V.ulm"))){

    if(!is.null(V$VBetween)){
      mBetween <- (abs(V$pBetween - V$qBetween) - 1) / 2
      nBetween <- (V$v - V$pBetween - 1) / 2
      FVBetween <- ((2 * nBetween + V$sBetween + 1) / (2 * mBetween + V$sBetween + 1)) * ((V$VBetween) / (V$sBetween - V$VBetween))

      termsPillaiBetween <- data.frame(df1 = V$sBetween * (2 * mBetween + V$sBetween + 1),
                                       df2 = V$sBetween * (2 * nBetween + V$sBetween + 1),
                                       VBetween = V$VBetween,
                                       FVBetween = FVBetween,
                                       p.value = 1 - pf(FVBetween, V$sBetween * (2 * mBetween + V$sBetween + 1),
                                                        V$sBetween * (2 * nBetween + V$sBetween + 1)))
      parmsPillaiBetween <- data.frame(p = V$pBetween, q = V$qBetween,
                                       s = V$sBetween,
                                       v = V$v, m = mBetween, n = nBetween)

      colnames(termsPillaiBetween) <- c("Num df", "Den df", "V", "F", "Pr(>F)")

      termsPillaiBetween <- cbind(termsPillaiBetween, mmemmuris::sigStars(termsPillaiBetween$`Pr(>F)`))

      colnames(termsPillaiBetween)[length(colnames(termsPillaiBetween))] <- ""

      if(any(class(V) %in% "V.ulm")){
        if(V$LR == TRUE){
          pillaiBetween <- list(pillai = termsPillaiBetween[(rownames(termsPillaiBetween) != "(Intercept)") & V$sBetween == 1, ],
                                parmsPillai = parmsPillaiBetween[(rownames(termsPillaiBetween) != "(Intercept)") & V$sBetween == 1, ],
                                approximation = "Pillai",
                                EB = V$EB, HB = V$HB,
                                Lbetween = V$Lbetween, L = V$L, type = V$type,
                                Vmatrix = V$EB / V$n)
        }else{
          pillaiBetween <- list(pillai = termsPillaiBetween[rownames(termsPillaiBetween) != "(Intercept)", ],
                                parmsPillai = parmsPillaiBetween[rownames(parmsPillaiBetween) != "(Intercept)", ],
                                approximation = "Pillai",
                                EB = V$EB, HB = V$HB,
                                Lbetween = V$Lbetween, L = V$L, type = V$type,
                                Vmatrix = V$EB / V$n)
        }
      }else{
        pillaiBetween <- list(pillai = termsPillaiBetween[rownames(termsPillaiBetween) != "(Intercept)", ],
                              parmsPillai = parmsPillaiBetween[rownames(parmsPillaiBetween) != "(Intercept)", ],
                              approximation = "Pillai",
                              EB = V$EB, HB = V$HB,
                              Lbetween = V$Lbetween, L = V$L, type = V$type,
                              Vmatrix = V$EB / V$n)
      }
    }else{
      termsPillaiBetween <- parmsPillaiBetween <- pillaiBetween <- NULL
    }

    mWithin <- (abs(V$pWithin - V$qWithin) - 1) / 2
    nWithin <- (V$v - V$pWithin - 1) / 2
    FVWithin <- ((2 * nWithin + V$sWithin + 1) / (2 * mWithin + V$sWithin + 1)) * ((V$VWithin) / (V$sWithin - V$VWithin))

    termsPillaiWithin <- data.frame(df1 = V$sWithin * (2 * mWithin + V$sWithin + 1),
                                    df2 = V$sWithin * (2 * nWithin + V$sWithin + 1),
                                    VWithin = V$VWithin,
                                    FVWithin = FVWithin,
                                    p.value = 1 - pf(FVWithin, V$sWithin * (2 * mWithin + V$sWithin + 1),
                                                     V$sWithin * (2 * nWithin + V$sWithin + 1)))
    parmsPillaiWithin <- data.frame(p = V$pWithin, q = V$qWithin,
                                    s = V$sWithin,
                                    v = V$v, m = mWithin, n = nWithin)

    colnames(termsPillaiWithin) <- c("Num df", "Den df", "V", "F", "Pr(>F)")

    termsPillaiWithin <- cbind(termsPillaiWithin, mmemmuris::sigStars(termsPillaiWithin$`Pr(>F)`))

    colnames(termsPillaiWithin)[length(colnames(termsPillaiWithin))] <- ""

    if(any(class(V) %in% "V.ulm")){
      if(V$LR == TRUE){
        pillaiWithin <- list(pillai = termsPillaiWithin[V$sWithin == 1, ],
                             parmsPillai = parmsPillaiWithin,
                             approximation = "Pillai",
                             E = V$E, H = V$H,
                             Lwithin = V$Lwithin, type = V$type)
      }else{
        pillaiWithin <- list(pillai = termsPillaiWithin,
                             parmsPillai = parmsPillaiWithin,
                             approximation = "Pillai",
                             E = V$E, H = V$H,
                             Lwithin = V$Lwithin, type = V$type)
      }
    }else{
      pillaiWithin <- list(pillai = termsPillaiWithin,
                           parmsPillai = parmsPillaiWithin,
                           approximation = "Pillai",
                           E = V$E, H = V$H,
                           Lwithin = V$Lwithin, type = V$type)
    }

    if(is.null(pillaiBetween))
      pillaiTests <- list(pillaiWithin = pillaiWithin)
    else
      pillaiTests <- list(pillaiBetween = pillaiBetween, pillaiWithin = pillaiWithin)
    if(any(class(V) %in% "V.mlm"))
      class(pillaiTests) <- c("list", "pillaiF", "mlm")
    else
      class(pillaiTests) <- c("list", "pillaiF", "ulm")
    return(pillaiTests)
  }else
    stop('Please provide an object of class "V.ulm" or "V.mlm".', call. = FALSE)
}

makePillaiTrace.ulm <- function(V, approximation = c("Muller", "Pillai")){
  approximation <- match.arg(approximation)

  if(approximation == "Muller"){
    mullerFApprox <- mullerF(V)

    withinTests <- mullerFApprox$pillaiWithin

    if(!is.null(mullerFApprox$pillaiBetween)){
      betweenTests <- mullerFApprox$pillaiBetween
      pillai <- list(betweenTests = betweenTests,
                     withinTests = withinTests,
                     LR = V$LR)
      class(pillai) <- c("list", "multiTest", "pillai", "ulm")
    }else{
      pillai <- list(withinTests = withinTests,
                     LR = V$LR)
      class(pillai) <- c("list", "multiTest", "pillai", "ulm")
    }
    return(pillai)
  }else if(approximation == "Pillai"){
    pillaiFApprox <- pillaiF(V)

    withinTests <- pillaiFApprox$pillaiWithin

    if(!is.null(pillaiFApprox$pillaiBetween)){
      betweenTests <- pillaiFApprox$pillaiBetween
      pillai <- list(betweenTests = betweenTests,
                     withinTests = withinTests,
                     LR = V$LR)
      class(pillai) <- c("list", "multiTest", "pillai", "ulm")
    }else{
      pillai <- list(withinTests = withinTests,
                     LR = V$LR)
      class(pillai) <- c("list", "multiTest", "pillai", "ulm")
    }
    return(pillai)
  }
}

mckeonF <- function(U){
  if(any(class(U) %in% c("U.mlm", "U.ulm"))){

    if(!is.null(U$UBetween)){
      nBetween <- (U$v - U$pBetween - 1) / 2
      bBetween <- ((U$pBetween + 2 * nBetween) * (U$qBetween + 2 * nBetween)) / (2 * (2 * nBetween + 1) * (nBetween - 1))
      CBetween <- (2 + (U$pBetween * U$qBetween + 2) / (bBetween - 1)) / (2 * nBetween)

      parmsHLBetween <- data.frame(p = U$pBetween, q = U$qBetween,
                                   s = U$sBetween, v = U$v,
                                   n = nBetween, b = bBetween, c = CBetween)

      FUBetween <- ((U$UBetween / CBetween) * ((4 + (U$pBetween * U$qBetween + 2) / (bBetween - 1)))) / (U$pBetween * U$qBetween)
      dfNumBetween <- U$pBetween * U$qBetween
      dfDenBetween <- 4 + (U$pBetween * U$qBetween + 2) / (bBetween - 1)

      termsHLBetween <- data.frame(df1 = dfNumBetween, df2 = dfDenBetween, U = U$UBetween, F = FUBetween,
                                   p.value = 1 - pf(FUBetween, dfNumBetween, dfDenBetween))

      colnames(termsHLBetween) <- c("Num df", "Den df", "U", "F", "Pr(>F)")
      termsHLBetween <- cbind(termsHLBetween, mmemmuris::sigStars(termsHLBetween$`Pr(>F)`))
      colnames(termsHLBetween)[length(colnames(termsHLBetween))] <- ""

      hltBetween <- list(hlt = termsHLBetween[rownames(termsHLBetween) != "(Intercept)", ],
                         parmsHL = parmsHLBetween[rownames(parmsHLBetween) != "(Intercept)", ],
                         approximation = "McKeon",
                         EB = U$EB, HB = U$HB,
                         Lbetween = U$Lbetween, L = U$L, type = U$type,
                         Vmatrix = U$EB / U$n)
    }else{
      termsHLBetween <- parmsHLBetween <- hltBetween <- NULL
    }

    nWithin <- (U$v - U$pWithin - 1) / 2
    bWithin <- ((U$pWithin + 2 * nWithin) * (U$qWithin + 2 * nWithin)) / (2 * (2 * nWithin + 1) * (nWithin - 1))
    CWithin <- (2 + (U$pWithin * U$qWithin + 2) / (bWithin - 1)) / (2 * nWithin)

    parmsHLWithin <- data.frame(p = U$pWithin, q = U$qWithin,
                                s = U$sWithin, v = U$v,
                                n = nWithin, b = bWithin, c = CWithin)

    FUWithin <- ((U$UWithin / CWithin) * ((4 + (U$pWithin * U$qWithin + 2) / (bWithin - 1)))) / (U$pWithin * U$qWithin)
    dfNumWithin <- U$pWithin * U$qWithin
    dfDenWithin <- 4 + (U$pWithin * U$qWithin + 2) / (bWithin - 1)

    termsHLWithin <- data.frame(df1 = dfNumWithin, df2 = dfDenWithin, U = U$UWithin, F = FUWithin,
                                p.value = 1 - pf(FUWithin, dfNumWithin, dfDenWithin))

    colnames(termsHLWithin) <- c("Num df", "Den df", "U", "F", "Pr(>F)")
    termsHLWithin <- cbind(termsHLWithin, mmemmuris::sigStars(termsHLWithin$`Pr(>F)`))
    colnames(termsHLWithin)[length(colnames(termsHLWithin))] <- ""

    hltWithin <- list(hlt = termsHLWithin[rownames(termsHLWithin) != "(Intercept)", ],
                      parmsHL = parmsHLWithin[rownames(parmsHLWithin) != "(Intercept)", ],
                      approximation = "McKeon",
                      E = U$E, H = U$H,
                      Lwithin = U$Lwithin, type = U$type)

    if(is.null(hltBetween))
      hltTests <- list(hltWithin = hltWithin)
    else
      hltTests <- list(hltBetween = hltBetween, hltWithin = hltWithin)
    if(any(class(U) %in% "U.mlm"))
      class(hltTests) <- c("list", "mckeonF", "mlm")
    else
      class(hltTests) <- c("list", "mckeonF", "ulm")
    return(hltTests)
  }else
    stop('Please provide an object of class "U.ulm" or "U.mlm".', call. = FALSE)
}

psF <- function(U){
  if(any(class(U) %in% c("U.mlm", "U.ulm"))){

    if(!is.null(U$UBetween)){
      mBetween <- (abs(U$pBetween - U$qBetween) - 1) / 2
      nBetween <- (U$v - U$pBetween - 1) / 2

      parmsHLBetween <- data.frame(p = U$pBetween, q = U$qBetween,
                                   s = U$sBetween, v = U$v, m = mBetween,
                                   n = nBetween)

      FUBetween <- (2 * (U$sBetween * nBetween + 1) * U$UBetween) / (U$sBetween ^ 2 * (2 * mBetween + U$sBetween + 1))
      dfNumBetween <- U$sBetween * (2 * mBetween + U$sBetween + 1)
      dfDenBetween <- 2 * (U$sBetween * nBetween + 1)

      termsHLBetween <- data.frame(df1 = dfNumBetween, df2 = dfDenBetween, U = U$UBetween, F = FUBetween,
                                   p.value = 1 - pf(FUBetween, dfNumBetween, dfDenBetween))

      colnames(termsHLBetween) <- c("Num df", "Den df", "U", "F", "Pr(>F)")
      termsHLBetween <- cbind(termsHLBetween, mmemmuris::sigStars(termsHLBetween$`Pr(>F)`))
      colnames(termsHLBetween)[length(colnames(termsHLBetween))] <- ""

      hltBetween <- list(hlt = termsHLBetween[rownames(termsHLBetween) != "(Intercept)", ],
                         parmsHL = parmsHLBetween[rownames(parmsHLBetween) != "(Intercept)", ],
                         approximation = "Pillai-Samson",
                         EB = U$EB, HB = U$HB,
                         Lbetween = U$Lbetween, L = U$L, type = U$type,
                         Vmatrix = U$EB / U$n)
    }else{
      termsHLBetween <- parmsHLBetween <- hltBetween <- NULL
    }

    mWithin <- (abs(U$pWithin - U$qWithin) - 1) / 2
    nWithin <- (U$v - U$pWithin - 1) / 2

    parmsHLWithin <- data.frame(p = U$pWithin, q = U$qWithin,
                                s = U$sWithin, v = U$v, m = mWithin,
                                n = nWithin)

    FUWithin <- (2 * (U$sWithin * nWithin + 1) * U$UWithin) / (U$sWithin ^ 2 * (2 * mWithin + U$sWithin + 1))
    dfNumWithin <- U$sWithin * (2 * mWithin + U$sWithin + 1)
    dfDenWithin <- 2 * (U$sWithin * nWithin + 1)

    termsHLWithin <- data.frame(df1 = dfNumWithin, df2 = dfDenWithin, U = U$UWithin, F = FUWithin,
                                p.value = 1 - pf(FUWithin, dfNumWithin, dfDenWithin))

    colnames(termsHLWithin) <- c("Num df", "Den df", "U", "F", "Pr(>F)")
    termsHLWithin <- cbind(termsHLWithin, mmemmuris::sigStars(termsHLWithin$`Pr(>F)`))
    colnames(termsHLWithin)[length(colnames(termsHLWithin))] <- ""

    hltWithin <- list(hlt = termsHLWithin[rownames(termsHLWithin) != "(Intercept)", ],
                      parmsHL = parmsHLWithin[rownames(parmsHLWithin) != "(Intercept)", ],
                      approximation = "Pillai-Samson",
                      E = U$E, H = U$H,
                      Lwithin = U$Lwithin, type = U$type)

    if(is.null(hltBetween))
      hltTests <- list(hltWithin = hltWithin)
    else
      hltTests <- list(hltBetween = hltBetween, hltWithin = hltWithin)
    if(any(class(U) %in% "U.mlm"))
      class(hltTests) <- c("list", "psF", "mlm")
    else
      class(hltTests) <- c("list", "psF", "ulm")
    return(hltTests)
  }else
    stop('Please provide an object of class "U.ulm" or "U.mlm".', call. = FALSE)
}

waldChi <- function(U){
  if(any(class(U) %in% c("U.mlm", "U.ulm"))){

    if(!is.null(U$UBetween)){
      parmsHLBetween <- data.frame(p = U$pBetween, q = U$qBetween, s = U$sBetween, v = U$v)
      termsHLBetween <- data.frame(df = U$pBetween * U$qBetween, U = U$UBetween,
                                   Chisq = U$X2Between,
                                   p.value = 1 - pchisq(U$X2Between, U$pBetween * U$qBetween))
      colnames(termsHLBetween) <- c("df", "U", "Chisq", "Pr(>Chisq)")
      termsHLBetween <- cbind(termsHLBetween, mmemmuris::sigStars(termsHLBetween$`Pr(>Chisq)`))
      colnames(termsHLBetween)[length(colnames(termsHLBetween))] <- ""

      hltBetween <- list(hlt = termsHLBetween[rownames(termsHLBetween) != "(Intercept)", ],
                         parmsHL = parmsHLBetween[rownames(parmsHLBetween) != "(Intercept)", ],
                         approximation = "Wald",
                         EB = U$EB, HB = U$HB,
                         Lbetween = U$Lbetween, L = U$L, type = U$type,
                         Vmatrix = U$EB / U$n)
    }else{
      termsHLBetween <- parmsHLBetween <- hltBetween <- NULL
    }

    parmsHLWithin <- data.frame(p = U$pWithin, q = U$qWithin, s = U$sWithin, v = U$v)
    termsHLWithin <- data.frame(df = U$pWithin * U$qWithin, U = U$UWithin,
                                Chisq = U$X2Within,
                                p.value = 1 - pchisq(U$X2Within, U$pWithin * U$qWithin))
    colnames(termsHLWithin) <- c("df", "U", "Chisq", "Pr(>Chisq)")
    termsHLWithin <- cbind(termsHLWithin, mmemmuris::sigStars(termsHLWithin$`Pr(>Chisq)`))
    colnames(termsHLWithin)[length(colnames(termsHLWithin))] <- ""

    hltWithin <- list(hlt = termsHLWithin[rownames(termsHLWithin) != "(Intercept)", ],
                      parmsHL = parmsHLWithin[rownames(parmsHLWithin) != "(Intercept)", ],
                      approximation = "Wald",
                      E = U$E, H = U$H,
                      Lwithin = U$Lwithin, type = U$type)

    if(is.null(hltBetween))
      hltTests <- list(hltWithin = hltWithin)
    else
      hltTests <- list(hltBetween = hltBetween, hltWithin = hltWithin)
    if(any(class(U) %in% "U.mlm"))
      class(hltTests) <- c("list", "waldChi", "mlm")
    else
      class(hltTests) <- c("list", "waldChi", "ulm")
    return(hltTests)
  }else
    stop('Please provide an object of class "U.ulm" or "U.mlm".', call. = FALSE)
}

modelCheck <- function(fMod){
  if(!(any(class(fMod) %in% c("gls", "lme"))))
    stop('The model must be a "gls" or "lme" object.', call. = FALSE)
}

getDataset <- function(fMod){
  if(any(class(fMod) %in% "lme")){
    if(is.null(fMod$data))
      stop('Please use the "data" argument in the model.  Make sure the "keep.data = FALSE" argument is not used.', call. = FALSE)
    dataset <- fMod$data
  }else{
    if(is.null(fMod$call$data))
      stop('Please use the "data" argument in the model.', call. = FALSE)
    dataset <- eval(fMod$call$data)
  }
  return(dataset)
}

#' Complete Cases
#'
#' @export
#'
#' @description
#' This function lets you know if you have any missing data in your dataset.
#'
#' @details
#' As introduced in Wright (1994) and (n.d.), the number of observations for
#' marginal and mixed-effects models is set equal to \deqn{n = N / r} where
#' \eqn{N} is the total number of observations of data in the long format and
#' \eqn{r} is the number of repeated measurements.  For multivariate linear models,
#' \eqn{N} is the total number of observations of data in the wide format and
#' \eqn{n} is the total number of complete observations of data in the wide format.
#'
#' @param fMod A model of class \code{\link[nlme]{gls}}, \code{\link[nlme]{lme}}, or \code{mlm}.
#'
#' @returns
#' \tabular{ll}{
#'    `complete`  \tab A logical telling the user if there are no missing values in the dataset. \cr
#'    `N`  \tab The total number of observations of data in the long or wide format. \cr
#'    `r`  \tab The total number of repeated measurements for data in the long format. \cr
#'    `n`  \tab The total number of observations of data used in the long or wide format. \cr
#' \tab \cr
#' }
#'
#' @references {Wright, S. P. (1994). Adjusted F Tests for Repeated Measures with
#' the MIXED Procedure. Knoxville: Statistics Department, University of Tennessee.
#' \url{https://support.sas.com/resources/papers/proceedings-archive/SUGI95/Sugi-95-195%20Wright.pdf}}
#' @references {Wright, S.P. (n.d.). Multivariate Analysis Using the MIXED Procedure. Sugi
#'   23: Multivariate analysis using the mixed procedure - SAS support.
#'   \url{https://support.sas.com/resources/papers/proceedings/proceedings/sugi23/Stats/p229.pdf}}
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#' completeData(glsMod)
#'
#' @seealso
#' {\code{\link[mmemmuris]{isReversible}}, \code{\link[mmemmuris]{BWInteracted}}}

completeData <- function(fMod){
  dataset <- mmemmuris:::getDataset(fMod)
  if(any(class(fMod) %in% c("gls", "lme"))){
      N <- nrow(nlme::getData(fMod))
     reps <- max(table(nlme::getGroups(fMod)))
     n <- N / reps
    complete <- (nrow(dataset) == N)
  }else{
    N <- nrow(dataset)
    n <- nrow(na.omit(dataset))
    reps <- ncol(mmemmuris::coefs(fMod))
    complete <- (N == n)
  }
  if(!complete){
    warning(
      paste0("Missing data is present. Number of experimental units set equal to ",
             n, "."),
      call. = FALSE
    )
  }
  if(any(class(fMod) %in% c("gls", "lme"))){
    return(list(complete = complete,
                N = N,
                r = reps,
                n = n))
  }else{
    return(list(complete = complete,
                N = N,
                n = n))
  }
}

#' Check If All Between-Subject Factors are Interacted with Within-Subject Factors
#'
#' @export
#'
#' @description
#' This function determines whether all between-subject terms are interacted
#' with the within-subject terms.
#'
#' @param fMod A model of class \code{\link[nlme]{gls}} or \code{\link[nlme]{lme}}.
#'
#' @returns
#' This function will return a logical (`TRUE`/`FALSE`).
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#' BWInteracted(glsMod)
#' glsMod <- gls(Distance ~ Sex + Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#' BWInteracted(glsMod)
#'
#' @seealso
#' {\code{\link[mmemmuris]{isReversible}}, \code{\link[mmemmuris]{completeData}}}

BWInteracted <- function(fMod){
  tt <- mmemmuris::termsType(fMod)
  orders <- attr(terms(fMod), "order")
  names(orders) <- attr(terms(fMod), "term.labels")
  withinOrders <- orders[(names(orders) %in% tt$within) & orders == 1]
  if(length(tt$between) == 0){
    if(any(class(fMod) == "gls"))
      fModFF <- update(fMod, model = as.formula(paste0(all.vars(terms(fMod))[1], paste0("~", names(withinOrders)))))
    else
      fModFF <- update(fMod, fixed = as.formula(paste0(all.vars(terms(fMod))[1], paste0("~", names(withinOrders)))))
  }else{
    betweenOrders <- orders[(names(orders) %in% tt$between)]
    if(any(class(fMod) == "gls"))
      fModFF <- update(fMod, model = as.formula(paste0(all.vars(terms(fMod))[1], paste0("~(", paste(names(betweenOrders), collapse = "+"), ")*", names(withinOrders)))))
    else
      fModFF <- update(fMod, fixed = as.formula(paste0(all.vars(terms(fMod))[1], paste0("~(", paste(names(betweenOrders), collapse = "+"), ")*", names(withinOrders)))))
  }
  termsFmodFF <- lapply(sapply(attr(terms(fModFF), "term.labels"), function(x) { strsplit(x, ":") }), function(x) { paste(sort(x), collapse = ":") })
  names(termsFmodFF) <- termsFmodFF
  termsFmod <- lapply(sapply(attr(terms(fMod), "term.labels"), function(x) { strsplit(x, ":") }), function(x) { paste(sort(x), collapse = ":") })
  names(termsFmod) <- termsFmod
  return((all(termsFmodFF[order(names(termsFmodFF))] %in% termsFmod[order(names(termsFmod))])) & (length(termsFmodFF[order(names(termsFmodFF))]) == length(termsFmod[order(names(termsFmod))])))
}

#' Check If Marginal and Mixed-Effects Models Correspond to a Multivariate Model
#'
#' @export
#'
#' @description
#' This function checks if a marginal (\code{\link[nlme]{gls}}) or mixed-effects
#' (\code{\link[nlme]{lme}}) model is reversible.  A marginal or mixed-effects
#' model is reversible if it can be equivalently specified as a multivariate linear model.
#'
#' @details
#' A reversible model needs complete cases (no missing data or all cases with
#' missing data omitted), all between-subject terms interacted with the
#' within-subject terms, and the marginal variance-covariance matrix of the
#' random/repeated effects for each subject must be the same unstructured
#' covariance structure.
#'
#' @param fMod A model of class \code{\link[nlme]{gls}} or \code{\link[nlme]{lme}}.
#'
#' @returns
#' This function returns a logical telling the user if the marginal or
#' mixed-effects model corresponds to a multivariate linear model.
#'
#' @references {\url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7009022/}}
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#'
#' # Reversible model
#' isReversible(glsMod)
#'
#' # Non-reversible model
#' isReversible(update(glsMod, model = Distance ~ Sex + Age))
#'
#' @seealso
#' {\code{\link[mmemmuris]{completeData}}, \code{\link[mmemmuris]{BWInteracted}}}

isReversible <- function(fMod){
  if(mmemmuris::covStruct(fMod) == "other")
    stop("Only unstructured and compound symmetry covariances are allowed.", call. = FALSE)
  return(mmemmuris::completeData(fMod)$complete & mmemmuris::BWInteracted(fMod))
}

determinantMatrix <- function(mat, tolerance = 1e-4){
  detMat <- det(mat)
  if(abs(detMat) < tolerance)
    return(0L)
  else
    return(detMat)
}

inverseMatrix <- function(mat, tolerance = 1e-4){
  invMat <- tryCatch(solve(mat, tol = tolerance), error = function(e) { e })
  if("error" %in% class(invMat))
    return(NULL)
  else
    return(invMat)
}

#' Epsilon Corrections for a Departure from Sphericity
#'
#' @export
#'
#' @description
#' This function will calculate Greenhouse-Geisser or Huynh-Feldt-Lecoutre
#' epsilons for univariate tests to correct within-subject effects for a
#' departure from sphericity for objects of class \code{\link[nlme]{gls}},
#' \code{\link[nlme]{lme}}, or \code{mlm}.
#'
#' @param Ematrices An object of class `Ematrix.ulm` or `Ematrix.mlm`.
#' @param method A character vector of which epsilon correction should be used.
#'
#' @returns
#' This function returns an epsilon value based on which correction was used.
#'
#' @details
#' The sum of squares cross products (SSCP) error matrices are
#' calculated as
#'
#' **E = M'** \eqn{(}**Y'Y**\eqn{-} **\eqn{\hat{B}}'** \eqn{(}**X'X**\eqn{)}**\eqn{\hat{B}}** \eqn{)}**M**
#'
#' where **M** is an identity matrix for between-subject effects and a
#' sum-to-zero contrast matrix for within-subject effects, **Y** is a combined
#' matrix of observations, **\eqn{\hat{B}}** is a matrix of combined coefficients from the
#' ordinary least squares models, and **X** is the design matrix of the
#' independent variables.
#'
#' For marginal and mixed-effects models, the SSCP **E** matrix is defined as
#'
#' **E =** \eqn{n}**M'VM**
#'
#' where \eqn{n} is the number of subjects, **M** is an identity matrix for
#' between-subject effects and a sum-to-zero contrast matrix for within-subject
#' effects, and **V** is the
#' marginal variance-covariance matrix where
#'
#' **V = ZGZ' + R** (**Z** is the random effects design matrix, **G** is the
#' corresponding (co)variances of the random effects matrix, and **R**
#' represents the (co)variances of the repeated effects matrix) for a
#' \code{\link[nlme]{gls}} or \code{\link[nlme]{lme}}
#' unstructured covariance model object (refit with maximum likelihood instead
#' of restricted maximum likelihood).
#'
#' The within-subject SSCP **E** matrix is used to calculate the
#' Greenhouse-Geisser and Huynh-Feldt-Lecoutre epsilons for univariate
#' adjustments for violations of sphericity.
#'
#' The Greenhouse-Geisser epsilon is calculated as
#'
#' \eqn{\epsilon_{GG} = ((\Sigma\lambda_i of} **E**\eqn{(}**M'M**\eqn{) ^ {-1}) ^ 2) / (r - 1)(\Sigma\lambda_i ^ 2 of} **E**\eqn{(}**M'M**\eqn{) ^ {-1})}
#'
#' The Huynh-Feldt-Lecoutre epsilon is calculated as
#'
#' \eqn{\epsilon_{HFL} = ((v + 1)(r - 1)\epsilon_{GG} - 2) / ((r - 1)(v - (r - 1)\epsilon_{GG}))}
#'
#' @references {\url{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.4/statug/statug_introreg_sect038.htm}}
#' @references {\url{https://www.lesahoffman.com/PSYC943/mv12psyc943_lecture13.pdf}}
#' @references {\url{https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.4/statug/statug_glm_details46.htm}}
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#'
#' # Greenhouse-Geisser epsilon
#' epsilon(Ematrix(glsMod))
#' # Huynh-Feldt-Lecoutre epsilon
#' epsilon(Ematrix(glsMod), method = "H")
#'
#' # Multivariate linear model
#' manovaMod <- manova(cbind(Distance8, Distance10, Distance12, Distance14) ~ Sex,
#'                     data = orthodontWide)
#'
#' # Greenhouse-Geisser epsilon
#' epsilon(Ematrix(manovaMod))
#' # Huynh-Feldt-Lecoutre epsilon
#' epsilon(Ematrix(manovaMod), method = "H")
#'
#' @seealso
#' {\code{\link[mmemmuris]{coefs}}, \code{\link[mmemmuris]{Vmatrix}}, \code{\link[mmemmuris]{Ematrix}}, \code{\link[mmemmuris]{sphericityTests}}}

epsilon <- function(Ematrices,
                    method = c("Greenhouse-Geisser", "Huynh-Feldt-Lecoutre")){
  method <- match.arg(method)
  if(any(class(Ematrices) %in% c("Ematrix.mlm", "Ematrix.ulm"))){
    if((Ematrices$r - 1) == 1L){
      warning("The lower-bound epsilon is 1.", call. = FALSE)
      return(1L)
    }
    if(is.null(mmemmuris:::inverseMatrix(Ematrices$E)))
      stop("The SSCP E matrix is singular.  Epsilon corrections are not available.", call. = FALSE)
    ggEpsilon <- (sum(diag(eigen(Ematrices$E %*% solve(t(Ematrices$M) %*% Ematrices$M))$values)) ^ 2) / ((Ematrices$r - 1) * sum(eigen(Ematrices$E %*% solve(t(Ematrices$M) %*% Ematrices$M))$values ^ 2))
    if(method == "Greenhouse-Geisser"){
      epsilon <- ggEpsilon
      if(ggEpsilon > 0.75)
        warning("The Huynh-Feldt-Lecoutre epsilon is preferred when the Greenhouse-Geisser epsilon is greater than 0.75.", call. = FALSE)
    }else{
      epsilon <- c(((Ematrices$n - rankMatrix(Ematrices$X) + 1) * (Ematrices$r - 1) * ggEpsilon - 2) / ((Ematrices$r - 1) * (Ematrices$n - rankMatrix(Ematrices$X) - (Ematrices$r - 1) * ggEpsilon)))
      if(ggEpsilon <= 0.75)
        warning("The Greenhouse-Geisser epsilon is preferred when it is less than or equal to 0.75.", call. = FALSE)
      if(any(min(epsilon, 1L) == 1L))
        warning("Huynh-Feldt-Lecoutre epsilon set to 1 for further calculations.", call. = FALSE)
    }
    return(epsilon)
  }else
    stop('Please provide an object of class "Ematrix.mlm" or "Ematrix.ulm".', call. = FALSE)
}

#' Likelihood Ratio Test with a Contrast Coefficients Matrix
#'
#' @export
#'
#' @description
#' This function performs a likelihood ratio test for a full and reduced model.
#'
#' @param fMod A model of class \code{\link[nlme]{gls}} or \code{\link[nlme]{lme}}.
#' @param L A contrast coefficients matrix.
#'
#' @returns
#' This function will return a `data.frame` for each term in the model.  Each
#' `data.frame` contains \tabular{ll}{
#'    `df` \tab The degrees of freedom used for the models. \cr
#'    \tab \cr
#'    `AIC` \tab The Akaike information criterions (AICs) for the models. \cr
#'    \tab \cr
#'    `BIC` \tab The Bayesian information criterions (BICs) for the models. \cr
#'    \tab \cr
#'    `logLik` \tab The log-likelihoods for the models. \cr
#'    \tab \cr
#'    `L.Ratio` \tab The likelihood ratio for the models. \cr
#'    \tab \cr
#'    `Pr(>Chisq)` \tab The p-value for the likelihood ratio test for the models. \cr
#' }
#'
#' @details
#' The maximum log-likelihood function for the full model is
#'
#' \eqn{logLik(M_F) = (-1 / 2)log(det(}**V_{F}** \eqn{)) - ((1 / 2)}**r'V_{F}** \eqn{^ {-1}}**r** \eqn{) - (N / 2)log(2\pi)}
#'
#' where **V_{F}** is the marginal variance-covariance matrix of the random/repeated
#' effects, **r** is the vector of model residuals, and \eqn{N} is the number of
#' observations of data in the long format.
#'
#' The contrast coefficients matrix (**L**) is calculated from the **LU**
#' decomposition of the crossproducts matrix (**X_{F}'X_{F}**).  To
#' project the design matrix of the reduced model (**X_{R}**) onto the
#' column space of the design matrix of the full model, the design matrix of the
#' full model needs to be multiplied by the orthogonal complement of **L'**.
#' It follows that the maximum log-likelihood function for the reduced model is
#'
#' \eqn{logLik(M_R) = (-1 / 2)log(det(}**V_{R}** \eqn{)) - ((1 / 2)}**r'V_{R}** \eqn{^ {-1}}**r** \eqn{) - (N / 2)log(2\pi)}
#'
#' where **V_{R}** is the marginal variance-covariance matrix of the random/repeated
#' effects, **r** is the vector of model residuals, and \eqn{N} is the number of
#' observations of data in the long format.
#'
#' The likelihood ratio (LR) statistic approximately follows a \eqn{\chi ^ 2}
#' distribution \deqn{LR = -2logLik(M_R) - {-2logLik(M_F)}} with
#' \eqn{df_{M_R} - df_{M_F}} degrees of freedom.
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#'
#' # Likelihood ratio chi-square statistic
#' L <- mmemmuris:::Lmatrix.ulm(glsMod)
#' L$Sex <- rbind(L$Sex, L$`Sex:Age`)
#' LRT(glsMod, L)
#'
#' # Likelihood ratio chi-square statistic
#' wilksLambda(glsMod, approximation = "LR")
#'
#' @seealso
#' {\code{\link[mmemmuris]{REMLtoML}}, \code{\link[mmemmuris]{wilksLambda}}}

LRT <- function(fMod, L){
  fModML <- mmemmuris::REMLtoML(fMod)
  X <- list(model.matrix(fModML, data = na.omit(mmemmuris:::getDataset(fModML))))
  rModX <- Map(function(x, y) { mmemmuris:::colSpaceX2inColSpaceX1(x, y) }, X, L)
  formulae <- lapply(rModX, function(x){
    as.formula(paste(rownames(attr(terms(fMod), "factors"))[1], "~ 0 +", paste(colnames(x), collapse = "+")))
  })
  dat <- lapply(rModX, function(x){
    na.omit(cbind(x, na.omit(mmemmuris:::getDataset(fMod))))
  })
  rModML <- Map(function(x, y) { update(fModML, eval(x), data = y) }, formulae, dat)
  LR <- lapply(rModML, function(x){
    LR <- anova(x, fModML, test = TRUE)
    rownames(LR) <- c("rModML", "fModML")
    LR
  })
  names(LR) <- names(L)
  return(LR)
}

#' Wald-Type F-Test with a Contrast Coefficients Matrix
#'
#' @export
#'
#' @description
#' This function calculates Wald-type F-statistics and Wald chi-square
#' statistics for a model.
#'
#' @param fMod A model of class \code{\link[nlme]{gls}} or \code{\link[nlme]{lme}}.
#' @param L A contrast coefficients matrix.
#'
#' @returns \tabular{ll}{
#'    `F` \tab A vector of Wald-type F-statistics. \cr
#'    \tab \cr
#'    `NumDF` \tab A vector of numerator degrees of freedom. \cr
#'    \tab \cr
#'    `X2` \tab A vector of Wald chi-square statistics.
#'    \cr \tab \cr
#' }
#'
#' @details
#' For marginal and mixed-effects models, we start with a Wald-type quadratic
#' form
#'
#' \eqn{F = (}**\eqn{\hat{\beta}}'L**\eqn{(}**L'** \eqn{(}**X'V**\eqn{ ^ {-1}}**X**\eqn{) ^ {-}}**L**\eqn{) ^ {-1}}**L'\eqn{\hat{\beta}}** \eqn{) / rank(}**L'** \eqn{(}**X'V**\eqn{ ^ {-1}}**X**\eqn{) ^ {-}}**L**\eqn{)}
#'
#' where **\eqn{\hat{\beta}}** is a vector of fixed effects coefficients, **X** is the fixed effects design matrix, and
#' **V** is the marginal variance-covariance matrix
#' where **V = ZGZ' + R** (**Z** is the random effects design matrix, **G** is
#' the corresponding (co)variances of the random effects matrix, and **R**
#' represents the (co)variances of the repeated effects matrix).  The contrast coefficients matrix (**L**) is calculated from the **LU**
#' decomposition of the crossproducts matrix (**X'X**).
#'
#' It follows that the Wald-\eqn{\chi ^ 2} statistic is the numerator degrees of
#' freedom (\eqn{rank(}**L'** \eqn{(}**X'V**\eqn{ ^ {-1}}**X**\eqn{) ^ {-}}**L**\eqn{)})
#' multiplied by the F-statistic.
#'
#' @examples
#' # Marginal model
#' # Unstructured covariance
#' glsMod <- gls(Distance ~ Sex * Age,
#'               correlation =
#'                 corSymm(form = ~ as.numeric(Age) | Subject),
#'               weights = varIdent(form = ~ 1 | Age),
#'               na.action = na.omit,
#'               data = orthodontLong)
#'
#' # Wald chi-square statistic
#' L <- mmemmuris:::Lmatrix.ulm(glsMod)
#' L$Sex <- rbind(L$Sex, L$`Sex:Age`)
#' waldF(glsMod, L)
#'
#' # Wald chi-square statistic
#' hlTrace(glsMod, approximation = "W")
#'
#' # Wald chi-square statistic
#' Ftests <- anova(glsMod)[c("Age", "Sex:Age"), ]
#' Ftests$Chisq <- Ftests$numDF * Ftests$`F-value`
#' Ftests$`Pr(>Chisq)` <- 1 - pchisq(Ftests$Chisq, Ftests$numDF)
#' Ftests
#'
#' @seealso
#' \code{\link[mmemmuris]{hlTrace}}

waldF <- function(fMod, L){
  if(any(class(fMod) %in% c("gls", "lme"))){
    if(any(class(fMod) %in% "lme"))
      # Don't like this workaround but seems to work
      fMod$fixDF$X[names(fMod$fixDF$X)] <- mmemmuris::ddfBW(fMod)$wDenDF
    wF <- lapply(L, function(x){ anova(fMod, L = x) })
    NumDF <- sapply(wF, function(x){ x$numDF })
    F <- sapply(wF, function(x){ x$`F-value` })
    X2 <- NumDF * F
    return(list(F = F, NumDF = NumDF, X2 = X2))
  }
  else
    stop('Please provide a model of class "gls" or "lme".', call. = FALSE)
}

LRT.cs <- function(fMod, rMod){
  datFull <- na.omit(mmemmuris:::getDataset(fMod))
  datRed <- na.omit(mmemmuris:::getDataset(rMod))
  VFull <- (fMod$df.residual / nrow(datFull)) * mmemmuris::Vmatrix(fMod)$csV
  VFull <- as.matrix(bdiag(replicate(nrow(datFull), VFull, simplify = FALSE)))
  rFull <- cbind(c(t(residuals(fMod))))
  NFull <- nrow(datFull) * length(colnames(coef(fMod)))
  csFull <- (-1 / 2) * log(det(VFull)) - ((1 / 2) * t(rFull) %*% solve(VFull) %*% rFull) - (NFull / 2) * log(2 * pi)
  VRed <- (rMod$df.residual / nrow(datRed)) * mmemmuris::Vmatrix(rMod)$csV
  VRed <- as.matrix(bdiag(replicate(nrow(datRed), VRed, simplify = FALSE)))
  rRed <- cbind(c(t(residuals(rMod))))
  NRed <- nrow(datRed) * length(colnames(coef(rMod)))
  csRed <- (-1 / 2) * log(det(VRed)) - ((1 / 2) * t(rRed) %*% solve(VRed) %*% rRed) - (NRed / 2) * log(2 * pi)
  dfModel <- c(length(coef(rMod)) + 2, length(coef(fMod)) + 2)
  X2 <- -2 * csRed - (-2 * csFull)
  lik <- c(csRed, csFull)
  N <- c(NRed, NFull)
  nestedModels <- data.frame(df = dfModel, AIC = 2 * dfModel - 2 * lik,
                             BIC = dfModel * log(N) - 2 * lik, logLik = lik,
                             L.Ratio = c(NA, X2),
                             p.value = c(NA, 1 - pchisq(X2, diff(dfModel))))
  rownames(nestedModels) <- c("rMod", "fMod")
  return(nestedModels)
}
