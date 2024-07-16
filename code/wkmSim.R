#' @title Wasserstein-Kaplan-Meier Survival Regression
#' @description Wrapper function for simulations.
#' @param y A matrix of quantile functions. Each row gives a quantile function.
#' @param x An matrix of predictors. Each row gives a \eqn{p}-dimensional 
#' predictor for the corresponding quantile function in \code{y}.
#' @param xOut A matrix of output predictor levels.
#' @param optns a list of options control parameters specified by
#' \code{list(name = value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{lower}{The smallest possible survival time. Default is 0.}
#' \item{upper}{The largest possible survival time. Default is `Inf`.}
#' }
#' @return A list containing the following fields:
#' \item{qf}{A matrix of quantile functions corresponding to \code{xOut}. 
#' Each row contains the quantile function of the survival time for a predctor level.}
#' \item{qfsupp}{A numeric vector representing the domain grid of quantile 
#' functions in `qf`.}
#' \item{y}{The quantile functions used.}
#' \item{x}{The predictors used.}
#' \item{xOut}{The output predictor levels used.}
#' \item{optns}{The control options used.}
#' @references
#' \itemize{
#' \item \cite{Zhou, Y. and Müller, H.-G. (2024+). Wasserstein-Kaplan-Meier Survival Regression.}
#' }

wkm <- function(y = NULL,
                x = NULL,
                xOut = NULL,
                optns = list()) {
  if(is.null(optns$lower)) {
    optns$lower <- 0
  }
  n <- nrow(x)
  nOut <- nrow(xOut)
  M <- ncol(y)
  
  # initialization of OSQP solver
  A <- cbind(diag(M), rep(0, M)) + cbind(rep(0, M), -diag(M))
  if (!is.null(optns$upper)) {
    # if optns$upper is not NULL
    l <- c(optns$lower, rep(0, M - 1), -optns$upper)
  } else {
    # if optns$upper is NULL
    A <- A[, -ncol(A)]
    l <- c(optns$lower, rep(0, M - 1))
  }
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
  } # for p = 1
  
  qf <- matrix(nrow = nOut, ncol = M)
  for (i in 1:nOut) {
    w <- apply(wc, 1, function(wci) {
      1 + t(wci) %*% (xOut[i, ] - xMean)
    })
    qNew <- apply(y, 2, weighted.mean, w) # M
    if (any(w < 0)) {
      # if negative weights exist
      model$Update(q = -qNew)
      qNew <- sort(model$Solve()$x)
    }
    if (!is.null(optns$upper)) {
      qNew <- pmin(qNew, optns$upper)
    }
    qNew <- pmax(qNew, optns$lower)
    qf[i, ] <- qNew
  }
  qfsupp <- 1:M / M
  res <-
    list(
      qf = qf,
      qfsupp = qfsupp,
      y = y,
      x = x,
      xOut = xOut,
      optns = optns
    )
  res
}