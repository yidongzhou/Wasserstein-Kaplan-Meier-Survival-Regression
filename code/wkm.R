#' @title Wasserstein-Kaplan-Meier Survival Regression
#' @description Performs Wasserstein-Kaplan-Meier survival regression for heterogeneous populations.
#' @param df A data frame containing the following columns:
#'   \describe{
#'     \item{time}{Survival time.}
#'     \item{censor}{Censoring indicator, where 0 indicates censored observations.}
#'     \item{...}{Additional covariate columns.}
#'   }
#' @param optns A list of control options specified as \code{list(name = value)}. 
#' See `Details'.
#' @details
#' Available control options:
#'   \describe{
#'     \item{lower}{The smallest possible survival time. Default is 0.}
#'     \item{upper}{The largest possible survival time. Default is \code{Inf}.}
#'   }
#' @return A \code{wkm} object, which is a list containing the following components:
#'   \describe{
#'     \item{subgroup}{A matrix of covariate values for each subgroup.}
#'     \item{qf}{A matrix of quantile functions corresponding to \code{subgroup}. 
#'     Each row contains the quantile function of the survival time for a subgroup.}
#'     \item{qfsupp}{A numeric vector representing the domain grid of quantile 
#'     functions in \code{qf}.}
#'     \item{df}{The input data frame.}
#'     \item{optns}{The control options used.}
#'   }
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
  if(is.null(optns$lower)) {
    optns$lower <- 0
  }
  cov <- setdiff(names(df), c('time', 'censor'))
  x0 <- as.matrix(expand.grid(lapply(cov, function(covi) 1:length(unique(df[, covi])) - 1)))
  colnames(x0) <- cov
  dfl <- split(df, df[, cov])
  idx <- sapply(dfl, function(dfli) sum(dfli$censor)) > 0
  dfl <- dfl[idx]# exclude groups with no events
  x <- x0[idx, ]
  n0 <- nrow(x0)
  n <- nrow(x) # number of observations
  
  # obtain Kaplan-Meier quantile functions
  M <- 5000
  y <- matrix(nrow = n, ncol = M)
  for(i in 1:n) {
    sc <- survfit(Surv(time, censor) ~ 1, data = dfl[[i]])
    obs <- cbind(sc$time, sc$surv)
    wt <- diff(c(0, 1 - obs[, 2]))
    sp <- c(obs[wt != 0, 1], max(dfl[[i]]$time))# Efron
    wt <- wt[wt != 0]
    yi <- rep(sp, M * c(wt, 1 - sum(wt)))
    if(length(yi) < M) {
      yi <- c(yi, sample(sp, M - length(yi), prob = c(wt, 1 - sum(wt))))
    }
    y[i, ] <- sort(yi)
  }
  
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
  
  qf <- matrix(nrow = n0, ncol = M)
  for (i in 1:n0) {
    w <- apply(wc, 1, function(wci) {
      1 + t(wci) %*% (x0[i, ] - xMean)
    })
    qNew <- apply(y, 2, weighted.mean, w) # M
    if (any(w < 0)) {
      # if negative weights exist, project
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
  
  res <- list(
    subgroup = x0,
    qf = qf,
    qfsupp = qfsupp,
    df = df,
    optns = optns
  )
  
  class(res) <- "wkm"
  res
}