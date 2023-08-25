
# code/TAGC/TAGC/TAGC_fisher_CI.r

CorCI_fz<-function(rho, n, sigma0, conf.level = 0.95, alternative = c("two.sided", "less", "greater")){
  if (n < 3L)
    stop("not enough finite observations")
  if ( (!missing(conf.level) && (length(conf.level) != 1) || !is.finite(conf.level) || conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  alternative <- match.arg(alternative)
  if (identical(rho, 1))
    ci <- c(1, 1)
  else {
    z <- DescTools::FisherZ(rho)
    sigma <-sigma0
    ci <- switch(alternative,
                 less = c(-Inf, z + sigma * qnorm(conf.level)),
                 greater = c(z - sigma * qnorm(conf.level), Inf),
                 two.sided = z + c(-1, 1) * sigma * qnorm((1 + conf.level)/2))
    ci <- DescTools::FisherZInv(ci)
  }
  return(c(cor = rho, lwr.ci = ci[1], upr.ci = ci[2]))
}

