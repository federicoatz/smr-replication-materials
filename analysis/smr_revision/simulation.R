# Run from the project root. Outputs are isolated from the original analysis.
suppressPackageStartupMessages(library(pwr))
set.seed(20260916)
out <- 'analysis/smr_revision'
n_for_r <- function(r) {
  if (pwr.r.test(n=4,r=abs(r),sig.level=.05)$power >= .8) return(4)
  tryCatch(ceiling(pwr.r.test(r=abs(r),sig.level=.05,power=.8)$n), error=function(e) NA_real_)
}
settings <- expand.grid(rho=c(.1,.3), median_N=c(100,500), selected=c(FALSE,TRUE))
results <- list(); diagnostics <- list(); k <- 0
for (s in seq_len(nrow(settings))) {
  cfg <- settings[s,]
  for (rep in 1:100) {
    N <- pmax(10,round(exp(rnorm(200,log(cfg$median_N),.7))))
    r <- tanh(rnorm(200,atanh(cfg$rho),1/sqrt(N-3)))
    t <- abs(r)*sqrt((N-2)/(1-r^2))
    keep <- t > if (cfg$selected) 1.96 else qnorm(.8)
    N <- N[keep]; r <- r[keep]; t <- t[keep]
    naive <- vapply(r,n_for_r,numeric(1))
    safe <- vapply(r*(1-qnorm(.8)/t),n_for_r,numeric(1))
    finite <- is.finite(naive) & is.finite(safe)
    diagnostics[[length(diagnostics)+1]] <- data.frame(cfg,rep=rep,retained=length(N),solver_excluded=sum(!finite))
    N <- N[finite]; naive <- naive[finite]; safe <- safe[finite]
    M <- safe/naive
    noise <- exp(rnorm(length(N),0,.15))
    rules <- list(point=naive,safeguard=safe,original_size=2*N,capped_safeguard=pmin(safe,500))
    for (rule in names(rules)) {
      # The same implementation noise is used across rules within each draw.
      observed <- pmax(4,round(if (rule == "capped_safeguard") pmin(safe*noise,500) else rules[[rule]]*noise))
      d <- data.frame(y=log(observed),x=log(naive),m=log(M),o=log(N))
      a <- lm(y~x+m,d); b <- lm(y~x+m+o,d)
      if (rule %in% c("point","safeguard")) {
        control <- d; control$y <- log(rules[[rule]])
        expected <- if (rule == "point") 0 else 1
        stopifnot(abs(coef(lm(y~x+m,control))["m"]-expected)<1e-7,
                  abs(coef(lm(y~x+m+o,control))["m"]-expected)<1e-7)
      }
      k <- k+1
      results[[k]] <- data.frame(cfg,rep=rep,rule=rule,n=length(N),
        portfolio=unname(coef(a)['m']),conditioned=unname(coef(b)['m']),
        condition_number=kappa(model.matrix(b),exact=TRUE),attainment=mean(observed>=safe),median_ratio=median(observed/safe))
    }
  }
}
write.csv(do.call(rbind,diagnostics),file.path(out,"simulation_diagnostics.csv"),row.names=FALSE)
x <- do.call(rbind,results)
stopifnot(all(is.finite(x$portfolio)),all(is.finite(x$conditioned)))
write.csv(x,file.path(out,'simulation_draws.csv'),row.names=FALSE)
summary <- aggregate(cbind(portfolio,conditioned,attainment,median_ratio)~rho+median_N+selected+rule,x,mean)
summary$mcse_portfolio <- aggregate(portfolio~rho+median_N+selected+rule,x,function(z)sd(z)/sqrt(length(z)))$portfolio
summary$mcse_conditioned <- aggregate(conditioned~rho+median_N+selected+rule,x,function(z)sd(z)/sqrt(length(z)))$conditioned
write.csv(summary,file.path(out,'simulation_summary.csv'),row.names=FALSE)
print(aggregate(cbind(portfolio,conditioned,attainment)~rule,x,range))
# Paired expected slope differences check that target-driven rules are distinguishable
# in this controlled experiment; this is not a classifier for SCORE decisions.
a <- subset(x,rule=='point'); b <- subset(x,rule=='safeguard')
stopifnot(all(abs((b$portfolio-a$portfolio)-1)<.12),
          all(abs((b$conditioned-a$conditioned)-1)<.12))
cat('Paired target-rule slope checks passed.\n')
