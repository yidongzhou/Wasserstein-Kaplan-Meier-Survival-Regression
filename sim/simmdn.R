source("code/wkmSim.R")

# Parallel computing
library(doSNOW)
cl <- makeCluster(6)
registerDoSNOW(cl)
clusterEvalQ(cl, library(survival))
progress <- function(q) {
  if(q %% 10 == 0){
    cat(sprintf("%d runs are complete\n", q))
  }
}

Q <- 100
M <- 5000

p <- 0.5# success probability
nVec <- c(100, 200, 500)# sample size
# rhoVec <- c(0.05, 0.1, 0.5)# random variation
xOut <- data.frame(X = 0:1)

# generate random sample sizes N_i
set.seed(1)
N <- list()
for(i in 1:length(nVec)) {
  lambdan <- nVec[i]
  cn <- 0.5
  Ni <- NULL
  while (length(Ni) < nVec[i]) {
    Ni <- c(Ni, rpois(1, lambda = cn * lambdan))
    if (Ni[length(Ni)] < 1) {
      Ni <- Ni[-length(Ni)]
    }
  }
  N[[i]] <- Ni
}

####################
### Survival MDN ###
####################
start_time <- Sys.time()
set.seed(1)
yMean <- rbind(
  qexp(1:(M - 1) / M, rate = 2),
  sqrt(qexp(1:(M - 1) / M, rate = 2))
)
sew <- foreach(n = nVec, .combine = 'cbind') %:%
  foreach(icount(Q), .combine = 'c', .options.snow = list(progress = progress)) %dopar% {
    idx <- which(nVec == n)
    x <- matrix(rbinom(n, 1, p), ncol = 1)
    dt <- list()
    y <- matrix(nrow = n, ncol = M)
    for (i in 1:n) {
      lambda <- 2
      if(x[i, ]) {
        u <- sqrt(rexp(N[[idx]][i], rate = lambda))
      } else {
        u <- rexp(N[[idx]][i], rate = lambda)
      }
      cs <- runif(N[[idx]][i], min = 0, max = 2)
      while(all(u > cs)) {
        cs <- runif(N[[idx]][i], min = 0, max = 2)
      }
      dt[[i]] <- data.frame(
        time = pmin(u, cs),
        death = as.integer(u <= cs),
        X = rep(x[i], N[[idx]][i]))# X
      sc <- survfit(Surv(time, death) ~ 1, data = dt[[i]])
      obs <- cbind(sc$time, sc$surv)
      wt <- diff(c(0, 1 - obs[, 2]))
      sp <- c(obs[wt != 0, 1], max(dt[[i]]$time))# Efron
      wt <- wt[wt != 0]
      yi <- rep(sp, M * c(wt, 1 - sum(wt)))
      if (length(yi) < M) {
        yi <- c(yi, sample(sp, M - length(yi), prob = c(wt, 1 - sum(wt))))
      }
      y[i, ] <- sort(yi)
    }
    res <- wkm(y, x, as.matrix(xOut))
    
    sum((res$qf[, -M] - yMean) ^ 2) / (M - 1)# WKM
  }
end_time <- Sys.time()
end_time - start_time
save(sew, file = 'data/sew.RData')
stopCluster(cl)