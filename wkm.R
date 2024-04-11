#' @title Wasserstein-Kaplan-Meier Survival Regression
#' @description Wasserstein-Kaplan-Meier survival regression for heterogeneous populations.
#' @param df a data frame with columns \code{time}, \code{censor}, and other covariates. 
#' The column \code{time} is the survival time and \code{censor} is the censoring indicator 
#' where 0 indicates censored observations.
#' @param optns a list of options control parameters specified by
#' \code{list(name = value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{lower}{a scalar with the lower bound of the support of the measure. Default is \code{NULL}.}
#' \item{upper}{a scalar with the upper bound of the support of the measure. Default is \code{NULL}.}
#' }
#' @return A \code{wkm} object --- a list containing the following fields:
#' \item{qFit}{a matrix holding the quantile functions corresponding to \code{x}. 
#' Each row gives a quantile and the domain grid is given in \code{qFitSup}.}
#' \item{qFitSup}{a numeric vector giving the domain grid of \code{qFit}.}
#' \item{df}{the data frame used.}
#' \item{optns}{the control options used.}
#' @references
#' \itemize{
#' \item \cite{Zhou, Y. and Müller, H.-G. (2024+). Wasserstein-Kaplan-Meier Survival Regression.}
#' }
#' @export

wkm <- function(df = NULL, optns = list()) {
  library(survival)
  
  if(!is.data.frame(df)) {
    stop("df must be a data frame")
  }
  if(!all(c("time", "censor") %in% names(df))) {
    stop("df must contain both 'time' and 'censor'")
  }
  cov <- setdiff(names(df), c('time', 'censor'))
  x0 <- as.matrix(expand.grid(lapply(cov, function(covi) 1:length(unique(df[, covi])) - 1)))
  colnames(x0) <- cov
  dfl <- split(df, df[, cov])
  idx <- sapply(dfl, function(dfli) sum(dfli$censor)) > 0
  dfl <- dfl[idx]# exclude groups with no events
  x <- x0[idx, ]
  
  # obtain Kaplan-Meier quantile functions
  y <- list()
  for(i in 1:length(dfl)) {
    sc <- survfit(Surv(time, censor) ~ 1, data = dfl[[i]])
    obs <- cbind(sc$time, sc$surv)
    wt <- diff(c(0, 1 - obs[, 2]))
    sp <- c(obs[wt != 0, 1], max(dfl[[i]]$time))# Efron
    wt <- wt[wt != 0]
    y[[i]] <- rep(sp, M * c(wt, 1 - sum(wt)))
    if(length(y[[i]]) < M) {
      y[[i]] <- c(y[[i]], sample(sp, M - length(y[[i]]), prob = c(wt, 1 - sum(wt))))
    }
    y[[i]] <- sort(y[[i]])
  }
  
  n0 <- nrow(x0)
  n <- nrow(x) # number of observations
  
  M <- 5000
  yM <- matrix(unlist(y), ncol = M, byrow = TRUE) # n by M
  
  # initialization of OSQP solver
  A <- cbind(diag(M), rep(0, M)) + cbind(rep(0, M), -diag(M))
  if (!is.null(optns$upper) &
      !is.null(optns$lower)) {
    # if lower & upper are neither NULL
    l <- c(optns$lower, rep(0, M - 1), -optns$upper)
  } else if (!is.null(optns$upper)) {
    # if lower is NULL
    A <- A[, -1]
    l <- c(rep(0, M - 1), -optns$upper)
  } else if (!is.null(optns$lower)) {
    # if upper is NULL
    A <- A[, -ncol(A)]
    l <- c(optns$lower, rep(0, M - 1))
  } else {
    # if both lower and upper are NULL
    A <- A[, -c(1, ncol(A))]
    l <- rep(0, M - 1)
  }
  # P <- as(diag(M), "sparseMatrix")
  # A <- as(t(A), "sparseMatrix")
  P <- diag(M)
  A <- t(A)
  q <- rep(0, M)
  u <- rep(Inf, length(l))
  model <-
    osqp::osqp(
      P = P,
      q = q,
      A = A,
      l = l,
      u = u,
      osqp::osqpSettings(max_iter = 1e05, eps_abs = 1e-05, eps_rel = 1e-05, verbose = FALSE)
    )
  
  xMean <- colMeans(x)
  invVa <- solve(var(x) * (n - 1) / n)
  wc <-
    t(apply(x, 1, function(xi) {
      t(xi - xMean) %*% invVa
    })) # n by p
  if (nrow(wc) != n) {
    wc <- t(wc)
  } # for p=1
  
  qFit <- matrix(nrow = n0, ncol = M)
  for (i in 1:n0) {
    w <- apply(wc, 1, function(wci) {
      1 + t(wci) %*% (x0[i, ] - xMean)
    })
    qNew <- apply(yM, 2, weighted.mean, w) # M
    if (any(w < 0)) {
      # if negative weights exist, project
      model$Update(q = -qNew)
      qNew <- sort(model$Solve()$x)
    }
    if (!is.null(optns$upper)) {
      qNew <- pmin(qNew, optns$upper)
    }
    if (!is.null(optns$lower)) {
      qNew <- pmax(qNew, optns$lower)
    }
    qFit[i, ] <- qNew
  }
  qFitSup <- 1:M / M
  
  res <- list(
    qFit = qFit,
    qFitSup = qFitSup,
    df = df,
    optns = optns
  )
  
  class(res) <- "wkm"
  res
}