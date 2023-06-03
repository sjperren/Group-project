# simulate time series from state-space model
sim_Xt <- function(X0, n, gamma){

  # parameters
  phi <- exp(-1/gamma)
  sigma <- 1- exp(-2/gamma)

  # simulate time series of length n
  X <- c(X0)
  for(i in 1:n){
    Xi <- X[i]*phi + rnorm(1, 0, sigma)
    X <- c(X, Xi)
  }

  return(X)
}

#---------------------------------------------------------------------------------------------------------------------
# inference of GP under KF
gp.KF <- function(y, gamma, sigma.w, m0, sigma0){

  # compute the parameters
  phi <- exp(-1/gamma)
  sigma.v <- sqrt(1- exp(-2/gamma))

  # length of y
  n <- length(y)

  # record parameters m_{t|t}, sigma_{t|t}, m{t+1|t}, sigma_{t+1|t}
  tt.record <- data.frame(t = NULL, m = NULL, sigma = NULL)
  t1t.record <- data.frame(t = NULL, m = NULL, sigma = NULL)


  m1 <- sigma0^2*(y[1]-m0)/(sigma0^2+sigma.w^2)
  sigma1 <- sqrt(sigma0^2 - sigma0^4/(sigma0^2+sigma.w^2))

  tt.record1 <- data.frame(t = 1, m = m1, sigma = sigma1)
  tt.record <- rbind(tt.record, tt.record1)

  for(t in 1:(n-1)){
    # prediction step
    pred.param <- pred_step(m_tt = m1,
                            sigma_tt = sigma1,
                            phi = phi,
                            sigma_v = sigma.v)
    m0 <- pred.param$m
    sigma0 <- pred.param$sigma
    t1t.record.i <- data.frame(t = t, m = m0, sigma = sigma0)
    t1t.record <- rbind(t1t.record, t1t.record.i)

    # updating step
    update.param <-     update_step(m_t1t = m0,
                                    sigma_t1t = sigma0,
                                    phi = phi,
                                    sigma_w = sigma.w,
                                    new_y = y[t+1])
    m1 <- update.param$m
    sigma1 <- update.param$sigma
    tt.record.i <- data.frame(t = t+1, m = m1, sigma = sigma1)
    tt.record <- rbind(tt.record, tt.record.i)
  }
  # prediction step at t=n
  pred.param <- pred_step(m_tt = m1,
                          sigma_tt = sigma1,
                          phi = phi,
                          sigma_v = sigma.v)
  m0 <- pred.param$m
  sigma0 <- pred.param$sigma
  t1t.record.n <- data.frame(t = n, m = m0, sigma = sigma0)
  t1t.record <- rbind(t1t.record, t1t.record.n)

  return(list(tt=tt.record, t1t=t1t.record))
}

#---------------------------------------------------------------------------------------------------------------------
# smoothing under KF
smooth.KF <- function(n, tt.record, t1t.record, alpha){
  # alpha: (1-alpha) credible interval

  # initialize at t=n
  smooth.record <-tt.record %>% filter(t==n)

  for(s in (n-1):1){
    ss.param <- tt.record %>% filter(t==s)
    s1s.param <- t1t.record %>% filter(t==s)
    s1n.param <- smooth.record %>% filter(t==(s+1))

    Js <- phi*(ss.param$sigma/s1s.param$sigma)^2

    m.sn <- ss.param$m + Js*(s1n.param$m - s1s.param$m)
    sigma2.sn <- ss.param$sigma^2 + Js^2*(s1n.param$sigma^2 - s1s.param$sigma^2)

    smooth.record.s <- data.frame(t = s, m = m.sn, sigma = sqrt(sigma2.sn))
    smooth.record <- rbind(smooth.record, smooth.record.s)

  }

  # (1-alpha)% credible interval
  smooth.record$upper <- qnorm(1-alpha/2, mean=smooth.record$m, sd=smooth.record$sigma)
  smooth.record$lower <- qnorm(alpha/2, mean=smooth.record$m, sd=smooth.record$sigma)

  return(smooth.record)
}

#---------------------------------------------------------------------------------------------------------------------

loglikelihood.partial <- function(y, gamma, sigma.w, m0, sigma0){

  # inference of GP under KF
  KF.inference <- gp.KF(y, gamma, sigma.w, m0, sigma0)

  tt.record <- KF.inference$tt
  t1t.record <- KF.inference$t1t

  # initialise the data.frame
  partial.gamma.tt <- data.frame(t = NULL, m = NULL, sd2 = NULL)
  partial.sd2.tt <- data.frame(t = NULL, m = NULL, sd2 = NULL)

  partial.gamma.t1t <- data.frame(t = NULL, m = NULL, sd2 = NULL)
  partial.sd2.t1t <- data.frame(t = NULL, m = NULL, sd2 = NULL)

  partial.logp <- data.frame(k = NULL, partial.gamma = NULL, partial.sd2 = NULL)

  # compute base terms
  gmm <- 2*exp(-2/gamma)/gamma^2

  partial.gamma.21 <- data.frame(t = 1,
                                 m = -exp(-1/gamma)*y[1]/(1+sigma.w^2)/gamma^2,
                                 sd2 = 1 + gmm - gmm*(1-1/(1+sigma.w^2)))
  partial.gamma.t1t <- rbind(partial.gamma.t1t, partial.gamma.21)

  partial.sd2.21 <- data.frame(t = 1,
                               m = -exp(-1/gamma)*y[1]/(1+sigma.w^2)^2,
                               sd2 = exp(-2/gamma)/(1+sigma.w^2)^2)
  partial.sd2.t1t <- rbind(partial.sd2.t1t, partial.sd2.21)

  #loop of computing paritial derivative recursively
  n <- length(y)
  for(i in 1:(n-1)){
    partial.gamma.i <- partial.gamma.t1t %>% filter(t==i)
    partial.sd2.i <- partial.sd2.t1t %>% filter(t==i)
    t1t.record.i <- t1t.record %>% filter(t==i)
    tt.record.i <- tt.record %>% filter(t==(i+1))

    # some coefficients
    coef0 <- t1t.record.i$sigma^2 + sigma.w^2
    coef1 <- (y[i+1] - t1t.record.i$m)/coef0
    coef2 <- sigma.w^2/coef0

    # t = i+1
    # compute partial derivative of m_{t|t} and sigma.w_{t|t}

    # with respect to gamma
    partial.gamma.m.tt <- partial.gamma.i$m * coef2 + (y[i+1] - t1t.record.i$m)*partial.gamma.i$sd2*coef2/coef0
    partial.gamma.sd2.tt <- coef2^2 * partial.gamma.i$sd2

    # with respect to sigma.w^2
    partial.sd2.m.tt <- partial.sd2.i$m*coef2 + coef1/coef0*(partial.sd2.i$sd2*sigma.w^2 - t1t.record.i$sigma^2)

    partial.sd2.sd2.tt <- partial.sd2.i$sd2*coef2^2 + t1t.record.i$sigma^4/coef0^2

    # store the results
    partial.gamma.tt.i <- data.frame(t = i+1,
                                     m = partial.gamma.m.tt,
                                     sd2 = partial.gamma.sd2.tt)
    partial.gamma.tt <- rbind(partial.gamma.tt, partial.gamma.tt.i)
    partial.sd2.tt.i <- data.frame(t = i+1,
                                   m = partial.sd2.m.tt,
                                   sd2 = partial.sd2.sd2.tt)
    partial.sd2.tt <- rbind( partial.sd2.tt,  partial.sd2.tt.i)

    # compute partial derivative of m_{t+1|t} and sigma.w_{t+1|t}

    # with respect to gamma
    partial.gamma.m.t1t <- exp(-1/gamma)*tt.record.i$m/gamma^2 + exp(-1/gamma)*partial.gamma.m.tt
    partial.gamma.sd2.t1t <- -2*exp(-2/gamma)/gamma^2+ exp(-2/gamma)*partial.gamma.sd2.tt + 2/gamma^2*exp(-2/gamma)*tt.record.i$sigma^2

    # with respect to sigma.w^2
    partial.sd2.m.t1t <- exp(-1/gamma)*partial.sd2.m.tt
    partial.sd2.sd2.t1t <- exp(-2/gamma)*partial.sd2.sd2.tt

    partial.gamma.t1t.i <- data.frame(t = i+1,
                                      m = partial.gamma.m.t1t,
                                      sd2 = partial.gamma.sd2.t1t)
    partial.gamma.t1t <- rbind(partial.gamma.t1t, partial.gamma.t1t.i)
    partial.sd2.t1t.i <- data.frame(t = i+1,
                                    m = partial.sd2.m.t1t,
                                    sd2 = partial.sd2.sd2.t1t)
    partial.sd2.t1t <- rbind( partial.sd2.t1t,  partial.sd2.t1t.i)

    # log-likelihood term

    logp.gamma <-  coef1*partial.gamma.i$m -0.5*(1/coef0-coef1^2)*partial.gamma.i$sd2

    logp.sd2 <-  coef1*partial.sd2.i$m - 0.5*(1+partial.sd2.i$sd2)/coef0 + 0.5*coef1^2*(1+partial.sd2.i$sd2)

    partial.logp.i <- data.frame(k = i+1,
                                 partial.gamma = logp.gamma,
                                 partial.sd2 = logp.sd2)
    partial.logp <- rbind(partial.logp, partial.logp.i)
  }

  # compute sum of each log-likelihood term
  partial.gamma.logp <- sum(partial.logp$partial.gamma)
  partial.sd2.logp <- sum(partial.logp$partial.sd2) - 0.5/(1+sigma.w^2)


  return(list(gamma = partial.gamma.logp, sd2 = partial.sd2.logp))
}

#---------------------------------------------------------------------------------------------------------------------
# compute log-likelihood
compute_logp <- function(y, gamma, sigma.w, m0, sigma0){
  # compute the KF inference
  KF.result <- gp.KF(y, gamma, sigma.w, m0, sigma0)
  t1t.record <- KF.result$t1t

  # compute log-likelihood
  n <- length(y)
  y.mean <- t1t.record$m[-n]
  y.sd2 <- t1t.record$sigma[-n]^2+sigma.w^2
  value <- sum(log(dnorm(y[2:n], mean=y.mean, sd=sqrt(y.sd2)))) +
    log(dnorm(y[1], mean = m0, sd = sqrt(sigma0^2+sigma.w^2)))

  return(value)

}

# compute optimal hyper-parameters
optim_parm <- function(y, m0, sigma0){
  opt_param <- optim(par = c(1,1),
                     fn = function(parm) -1*compute_logp(y, parm[1],
                                                         parm[2], m0, sigma0))
  return(list(gamma = opt_param$par[1],
              sigma.w=opt_param$par[2]))

}

