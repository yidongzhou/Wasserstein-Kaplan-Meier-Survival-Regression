source("code/wkmSim.R")

# Parallel computing
library(doSNOW)
cl <- makeCluster(50)
registerDoSNOW(cl)
clusterEvalQ(cl, library(survival))
clusterEvalQ(cl, library(randomForestSRC))
progress <- function(q) {
  if(q %% 10 == 0){
    cat(sprintf("%d runs are complete\n", q))
  }
}

Q <- 100
M <- 5000

k <- 2# shape parameter
beta <- 1:5 / 100
p <- rep(0.5, 5)# success probability of each covariate
nVec <- c(50, 100, 200)# sample size
rhoVec <- c(0, 0.05, 0.1, 0.5)# random variation
ctVec <- c(2, 1)# censoring rate (2: 20%; 1: 50%)
xOut <- expand.grid(X1 = 0:1, X2 = 0:1, X3 = 0:1, X4 = 0:1, X5 = 0:1)

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

#####################
### Setting I WKM ###
#####################
sewr1 <- list()
serr1 <- list()
set.seed(1)
for(j in 1:length(rhoVec)) {
  rho <- rhoVec[j]
  sewr1[[j]] <- list()
  serr1[[j]] <- list()
  yMean <- t(apply(xOut, 1, function(xOuti)
    qweibull(1:(M - 1) / M, shape = k, scale = sum(xOuti * beta) + 0.1)))# conditional mean
  for(l in 1:length(ctVec)) {
    ct <- ctVec[l]
    se <- foreach(n = nVec, .combine = 'cbind') %:%
      foreach(icount(Q), .combine = 'c', .options.snow = list(progress = progress)) %dopar% {
        idx <- which(nVec == n)
        x <- sapply(p, function(pi) rbinom(n, 1, pi))# n by p
        # st <- Sys.time()
        dt <- list()
        y <- matrix(nrow = n, ncol = M)
        for (i in 1:n) {
          lambda <- ifelse(rho > 0, rgamma(1, shape = (sum(x[i, ] * beta) + 0.1)^2 / rho, 
                                           scale = rho / (sum(x[i, ] * beta) + 0.1)), 
                           sum(x[i, ] * beta) + 0.1)
          u <- rweibull(N[[idx]][i], shape = k, scale = lambda)
          cs <- rweibull(N[[idx]][i], shape = k, scale = ct * lambda)
          while(all(u > cs)) {
            cs <- rweibull(N[[idx]][i], shape = k, scale = ct * lambda)
          }
          dt[[i]] <- data.frame(
            time = pmin(u, cs),
            death = as.integer(u <= cs),
            t(replicate(N[[idx]][i], x[i, ])))# X1, ..., X5
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
        # et <- Sys.time()
        # et - st
        
        dt <- do.call('rbind', dt)
        obj <- rfsrc(Surv(time, death) ~ ., data = dt)
        # system.time(obj <- rfsrc(Surv(time, death) ~ ., data = dt))
        sc <- predict(obj, newdata = xOut)
        y0 <- list()
        wt0 <- list()
        for (i in 1:nrow(xOut)) {
          obs <- cbind(sc$time.interest, sc$survival[i, ])
          wt <- diff(c(0, 1 - obs[, 2]))
          y0i <- obs[wt != 0, 1]
          wt <- wt[wt != 0]
          if (sum(wt) >= 1) {
            if(length(wt) > 1) {
              y0i <- y0i[-length(y0i)]
              wt <- wt[-length(wt)]
            } else {# if a point mass
              y0i <- rep(y0i, M - 1)
              wt <- 1:(M - 1) / M
            }
          }
          y0[[i]] <- y0i
          wt0[[i]] <- wt
        }
        
        c(sum((res$qf[, -M] - yMean) ^ 2) / (M - 1),# WKM
          sum(sapply(1:nrow(xOut), function(i) {
            weighted.mean((y0[[i]] - qweibull(
              cumsum(wt0[[i]]),
              shape = k,
              scale = sum(xOut[i,] * beta) + 0.1
            )) ^ 2, w = wt0[[i]]) / sum(wt0[[i]])
          }))) / nrow(xOut)# RSF
      }
    sewr1[[j]][[l]] <- se[seq(1, 2 * Q, 2), ]# WKM
    serr1[[j]][[l]] <- se[seq(2, 2 * Q, 2), ]# RSF
    }
}
save(sewr1, serr1, file = 'data/ser1.RData')

######################
### Setting II CPH ###
######################
sewr2 <- list()
serr2 <- list()
set.seed(1)
for(j in 1:length(rhoVec)) {
  rho <- rhoVec[j]
  sewr2[[j]] <- list()
  serr2[[j]] <- list()
  yMean <- t(apply(xOut, 1, function(xOuti) 
    qweibull(1:(M - 1) / M, shape = k, scale = exp(-sum(xOuti * beta) / 2))))# conditional mean
  for(l in 1:length(ctVec)) {
    ct <- ctVec[l]
    se <- foreach(n = nVec, .combine = 'cbind') %:%
      foreach(icount(Q), .combine = 'c', .options.snow = list(progress = progress)) %dopar% {
        idx <- which(nVec == n)
        x <- sapply(p, function(pi) rbinom(n, 1, pi))# n by d
        dt <- list()
        y <- matrix(nrow = n, ncol = M)
        for (i in 1:n) {
          lambda <- ifelse(rho > 0, rgamma(1, shape = exp(-sum(x[i, ] * beta)) / rho, 
                                           scale = rho / exp(-sum(x[i, ] * beta) / 2)), 
                           exp(-sum(x[i, ] * beta) / 2))
          u <- rweibull(N[[idx]][i], shape = k, scale = lambda)
          cs <- rweibull(N[[idx]][i], shape = k, scale = ct * lambda)
          while(all(u > cs)) {
            cs <- rweibull(N[[idx]][i], shape = k, scale = ct * lambda)
          }
          dt[[i]] <- data.frame(
            time = pmin(u, cs),
            death = as.integer(u <= cs),
            t(replicate(N[[idx]][i], x[i, ])))# X1, ..., X5
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
        
        dt <- do.call('rbind', dt)
        obj <- rfsrc(Surv(time, death) ~ ., data = dt)
        sc <- predict(obj, newdata = xOut)
        y0 <- list()
        wt0 <- list()
        for (i in 1:nrow(xOut)) {
          obs <- cbind(sc$time.interest, sc$survival[i, ])
          wt <- diff(c(0, 1 - obs[, 2]))
          y0i <- obs[wt != 0, 1]
          wt <- wt[wt != 0]
          if (sum(wt) >= 1) {
            if(length(wt) > 1) {
              y0i <- y0i[-length(y0i)]
              wt <- wt[-length(wt)]
            } else {# if a point mass
              y0i <- rep(y0i, M - 1)
              wt <- 1:(M - 1) / M
            }
          }
          y0[[i]] <- y0i
          wt0[[i]] <- wt
        }
        
        c(sum((res$qf[, -M] - yMean) ^ 2) / (M - 1),# WKM
          sum(sapply(1:nrow(xOut), function(i) {
            weighted.mean((y0[[i]] - qweibull(
              cumsum(wt0[[i]]), shape = k, scale = exp(-sum(xOut[i,] * beta) / 2)
            )) ^ 2, w = wt0[[i]]) / sum(wt0[[i]])
          }))) / nrow(xOut)# RSF
      }
    sewr2[[j]][[l]] <- se[seq(1, 2 * Q, 2), ]# WKM
    serr2[[j]][[l]] <- se[seq(2, 2 * Q, 2), ]# RSF
  }
}
save(sewr2, serr2, file = 'data/ser2.RData')
stopCluster(cl)