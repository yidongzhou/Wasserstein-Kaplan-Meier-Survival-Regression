library(ggplot2)
Q <- 1000
nVec <- c(100, 200, 500)# sample size
rhoVec <- c(0.05, 0.1, 0.5)# random variation
ctVec <- c(2, 1)# censoring rate (2: 20%; 1: 50%)

# sim.R
load('data/se1.RData')
load('data/se2.RData')
for(j in 1:length(rhoVec)) {
  cat(paste(format(round((sapply(sew1[[j]], colMeans)), 4), scientific = FALSE), 
            collapse = '&'), '\n',
      paste(paste0('(', format(round((sapply(sew1[[j]], 
                                             function(sew1i) apply(sew1i, 2, sd))), 4), 
                               scientific = FALSE), ')'), collapse = '&'), '\n',
      paste(format(round((sapply(sec1[[j]], colMeans)), 4), scientific = FALSE), 
            collapse = '&'), '\n',
      paste(paste0('(', format(round((sapply(sec1[[j]], 
                                             function(sec1i) apply(sec1i, 2, sd))), 4), 
                               scientific = FALSE), ')'), collapse = '&'), '\n',
      paste(format(round((sapply(sew2[[j]], colMeans)), 4), scientific = FALSE), 
            collapse = '&'), '\n',
      paste(paste0('(', format(round((sapply(sew2[[j]], 
                                             function(sew2i) apply(sew2i, 2, sd))), 4), 
                               scientific = FALSE), ')'), collapse = '&'), '\n',
      paste(format(round((sapply(sec2[[j]], colMeans)), 4), scientific = FALSE), 
            collapse = '&'), '\n',
      paste(paste0('(', format(round((sapply(sec2[[j]], 
                                             function(sec2i) apply(sec2i, 2, sd))), 4), 
                               scientific = FALSE), ')'), collapse = '&'), '\n'
  )
}

# sim2.R
load('data/sem.RData')
load('data/sel.RData')
rhoVec <- c(0, 0.05, 0.1, 0.5)# random variation
for(j in 1:length(rhoVec)) {
  cat(paste(format(round((sapply(sewm[[j]], colMeans)), 4), scientific = FALSE), 
            collapse = '&'), '\n',
      paste(paste0('(', format(round((sapply(sewm[[j]], 
                                             function(sewmi) apply(sewmi, 2, sd))), 4), 
                               scientific = FALSE), ')'), collapse = '&'), '\n',
      paste(format(round((sapply(secm[[j]], colMeans)), 4), scientific = FALSE), 
            collapse = '&'), '\n',
      paste(paste0('(', format(round((sapply(secm[[j]], 
                                             function(secmi) apply(secmi, 2, sd))), 4), 
                               scientific = FALSE), ')'), collapse = '&'), '\n',
      paste(format(round((sapply(sewl[[j]], colMeans)), 4), scientific = FALSE), 
            collapse = '&'), '\n',
      paste(paste0('(', format(round((sapply(sewl[[j]], 
                                             function(sewli) apply(sewli, 2, sd))), 4), 
                               scientific = FALSE), ')'), collapse = '&'), '\n',
      paste(format(round((sapply(secl[[j]], colMeans)), 4), scientific = FALSE), 
            collapse = '&'), '\n',
      paste(paste0('(', format(round((sapply(secl[[j]], 
                                             function(secli) apply(secli, 2, sd))), 4), 
                               scientific = FALSE), ')'), collapse = '&'), '\n'
  )
}

# simb.R
load('data/seb.RData')
rhoVec <- c(0.05, 0.1, 0.5)# random variation
for(j in 1:length(rhoVec)) {
  cat(paste(format(round((sapply(sewb[[j]], colMeans)), 4), scientific = FALSE), 
            collapse = '&'), '\n',
      paste(paste0('(', format(round((sapply(sewb[[j]], 
                                             function(sewbi) apply(sewbi, 2, sd))), 4), 
                               scientific = FALSE), ')'), collapse = '&'), '\n',
      paste(format(round((sapply(secb[[j]], colMeans)), 4), scientific = FALSE), 
            collapse = '&'), '\n',
      paste(paste0('(', format(round((sapply(secb[[j]], 
                                             function(secbi) apply(secbi, 2, sd))), 4), 
                               scientific = FALSE), ')'), collapse = '&'), '\n'
  )
}

# simmdn.R
load('data/sew.RData')
round(colMeans(sew), 4)
round(apply(sew, 2, sd), 4)

# simrsf.R
load('data/ser1.RData')
load('data/ser2.RData')
rhoVec <- c(0, 0.05, 0.1, 0.5)# random variation
for(j in 1:length(rhoVec)) {
  cat(paste(format(round((sapply(sewr1[[j]], colMeans)), 4), scientific = FALSE), 
            collapse = '&'), '\n',
      paste(paste0('(', format(round((sapply(sewr1[[j]], 
                                             function(sewr1i) apply(sewr1i, 2, sd))), 4), 
                               scientific = FALSE), ')'), collapse = '&'), '\n',
      paste(format(round((sapply(serr1[[j]], colMeans)), 4), scientific = FALSE), 
            collapse = '&'), '\n',
      paste(paste0('(', format(round((sapply(serr1[[j]], 
                                             function(serr1i) apply(serr1i, 2, sd))), 4), 
                               scientific = FALSE), ')'), collapse = '&'), '\n',
      paste(format(round((sapply(sewr2[[j]], colMeans)), 4), scientific = FALSE), 
            collapse = '&'), '\n',
      paste(paste0('(', format(round((sapply(sewr2[[j]], 
                                             function(sewr2i) apply(sewr2i, 2, sd))), 4), 
                               scientific = FALSE), ')'), collapse = '&'), '\n',
      paste(format(round((sapply(serr2[[j]], colMeans)), 4), scientific = FALSE), 
            collapse = '&'), '\n',
      paste(paste0('(', format(round((sapply(serr2[[j]], 
                                             function(serr2i) apply(serr2i, 2, sd))), 4), 
                               scientific = FALSE), ')'), collapse = '&'), '\n'
  )
}

# simh.R
load('data/seh1.RData')
load('data/seh2.RData')
rhoVec <- c(0.05, 0.1, 0.5)# random variation
for(j in 1:length(rhoVec)) {
  cat(paste(format(round((sapply(sewh1[[j]], colMeans)), 4), scientific = FALSE), 
            collapse = '&'), '\n',
      paste(paste0('(', format(round((sapply(sewh1[[j]], 
                                             function(sewh1i) apply(sewh1i, 2, sd))), 4), 
                               scientific = FALSE), ')'), collapse = '&'), '\n',
      paste(format(round((sapply(sech1[[j]], colMeans)), 4), scientific = FALSE), 
            collapse = '&'), '\n',
      paste(paste0('(', format(round((sapply(sech1[[j]], 
                                             function(sech1i) apply(sech1i, 2, sd))), 4), 
                               scientific = FALSE), ')'), collapse = '&'), '\n',
      paste(format(round((sapply(sewh2[[j]], colMeans)), 4), scientific = FALSE), 
            collapse = '&'), '\n',
      paste(paste0('(', format(round((sapply(sewh2[[j]], 
                                             function(sewh2i) apply(sewh2i, 2, sd))), 4), 
                               scientific = FALSE), ')'), collapse = '&'), '\n',
      paste(format(round((sapply(sech2[[j]], colMeans)), 4), scientific = FALSE), 
            collapse = '&'), '\n',
      paste(paste0('(', format(round((sapply(sech2[[j]], 
                                             function(sech2i) apply(sech2i, 2, sd))), 4), 
                               scientific = FALSE), ')'), collapse = '&'), '\n'
  )
}

# boxplots
load('data/se1.RData')
load('data/se2.RData')
df1 <- data.frame(mspe = c(unlist(sew1), unlist(sec1)),
                 n = rep(rep(nVec, each = Q), length(rhoVec) * length(ctVec)),
                 rho = paste0('rho = ', rep(rhoVec, each = Q * length(ctVec) * length(nVec))),
                 method = factor(rep(c('WKM', 'CPH'), each = Q * length(rhoVec) * length(ctVec) * length(nVec)), levels = c('WKM', 'CPH')), 
                 censoring = rep(rep(paste0(c('20%', '50%'), ' censoring'), each = Q * length(nVec)), length(rhoVec)))
df2 <- data.frame(mspe = c(unlist(sew2), unlist(sec2)),
                  n = rep(rep(nVec, each = Q), length(rhoVec) * length(ctVec)),
                  rho = paste0('rho = ', rep(rhoVec, each = Q * length(ctVec) * length(nVec))),
                  method = factor(rep(c('WKM', 'CPH'), each = Q * length(rhoVec) * length(ctVec) * length(nVec)), levels = c('WKM', 'CPH')), 
                  censoring = rep(rep(paste0(c('20%', '50%'), ' censoring'), each = Q * length(nVec)), length(rhoVec)))
pdf(file = "latex/img/sim1.pdf", width = 13, height = 6.5)
ggplot(data = df1, aes(x = factor(n), y = mspe, color = method, linetype = method)) +
  geom_boxplot(outlier.alpha = 0.5) +
  ggh4x::facet_grid2(rows = vars(censoring), cols = vars(rho), scales = 'free_y', independent = 'y') +
  labs(x = 'n', y = 'MSPE') +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position = "top", legend.title = element_blank())
dev.off()
pdf(file = "latex/img/sim2.pdf", width = 13, height = 6.5)
ggplot(data = df2, aes(x = factor(n), y = mspe, color = method, linetype = method)) +
  geom_boxplot(outlier.alpha = 0.5) +
  ggh4x::facet_grid2(rows = vars(censoring), cols = vars(rho), scales = 'free_y', independent = 'y') +
  labs(x = 'n', y = 'MSPE') +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position = "top", legend.title = element_blank())
dev.off()
