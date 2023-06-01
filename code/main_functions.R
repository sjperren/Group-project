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

# prediction step
pred_step <- function(m.tt, sigma.tt, phi, sigma.v){
  
  # updated m_{t+1|t}
  m.t1t <- phi*m.tt
  
  # updated sigma_{t+1|t}
  sigma.t1t <- sqrt(sigma.v^2 + phi^2*sigma.tt^2)
  
  return(list(m=m.t1t, sigma=sigma.t1t))
}


# updating step
update_step <- function(m.t1t, sigma.t1t, phi, sigma.w, new.y){
  
  denominator <- sigma.w^2 + sigma.t1t^2
  
  # updated m_{t+1|t+1}
  m.t1t1 <- (sigma.w^2*m.t1t + sigma.t1t*new.y)/denominator
  
  # updated sigma_{t+1|t+1}
  sigma.t1t1 <- sqrt(sigma.w^2*sigma.t1t^2/denominator)
  
  return(list(m=m.t1t1, sigma=sigma.t1t1))
  
}



# KF inference
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
    pred.param <- pred_step(m1, sigma1, phi, sigma.v)
    m0 <- pred.param$m
    sigma0 <- pred.param$sigma
    t1t.record.i <- data.frame(t = t, m = m0, sigma = sigma0)
    t1t.record <- rbind(t1t.record, t1t.record.i)
    
    # updating step
    update.param <- update_step(m0, sigma0, phi, sigma.w, y[t+1])
    m1 <- update.param$m
    sigma1 <- update.param$sigma
    tt.record.i <- data.frame(t = t+1, m = m1, sigma = sigma1)
    tt.record <- rbind(tt.record, tt.record.i)
  }
  # prediction step at t=n
  pred.param <- pred_step(m1, sigma1, phi, sigma.v)
  m0 <- pred.param$m
  sigma0 <- pred.param$sigma
  t1t.record.n <- data.frame(t = n, m = m0, sigma = sigma0)
  t1t.record <- rbind(t1t.record, t1t.record.n)
  
  return(list(tt=tt.record, t1t=t1t.record))
}


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

# prediction under KF
pred.KF <- function(m.nn, sd.nn, gamma){
  # m.nn: numeric, smoothed m_{n|n}
  # sd.nn: numeric, smoothed sd_{n|n}
  
  # compute the parameters
  phi <- exp(-1/gamma)
  sigma.v <- sqrt(1- exp(-2/gamma))
  
  # predicted m
  pred.m <- phi*last.smooth$m
  
  # predicted sd^2
  pred.sd2 <- last.smooth$sigma^2*phi^2 + sigma.v^2
  
  # output 
  return(list(m = pred.m,
              sigma = sqrt(pred.sd2)))
}


loglikelihood.partial <- function(gamma, sigma.w){
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
  
  for(i in 1:(n-1)){
    partial.gamma.i <- partial.gamma.t1t %>% filter(t==i)
    partial.sd2.i <- partial.sd2.t1t %>% filter(t==i)
    t1t.record.i <- t1t.record %>% filter(t==i)
    tt.record.i <- tt.record %>% filter(t==i)
    
    
    partial.gamma.m.tt <- partial.gamma.i$m * sigma.w^2/(sigma.w^2 + t1t.record.i$sigma^2) +
      (y[i+1] - t1t.record.i$m)*partial.gamma.i$sd2*sigma.w^2/(sigma.w^2 + t1t.record.i$sigma^2)^2
    
    partial.gamma.sd2.tt <- sigma.w^4/(sigma.w^2 + t1t.record.i$sigma^2)^2 * partial.gamma.i$sd2
    
    partial.sd2.m.tt <- partial.sd2.i$m * sigma.w^2/(sigma.w^2 + t1t.record.i$sigma^2) + (y[i+1] - t1t.record.i$m)*(partial.sd2.i$sd2*sigma.w^2 - t1t.record.i$sigma^2)/(sigma.w^2+t1t.record.i$sigma^2)^2
    
    partial.sd2.sd2.tt <- partial.sd2.i$sd2*sigma.w^2/(sigma.w^2+t1t.record.i$sigma^2) - t1t.record.i$sigma^2*(partial.sd2.i$sd2*sigma.w^2 - t1t.record.i$sigma^2)/(sigma.w^2+t1t.record.i$sigma^2)^2
    
    
    partial.gamma.tt.i <- data.frame(t = i+1, 
                                     m = partial.gamma.m.tt, 
                                     sd2 = partial.gamma.sd2.tt)
    partial.gamma.tt <- rbind(partial.gamma.tt, partial.gamma.tt.i)
    partial.sd2.tt.i <- data.frame(t = i+1, 
                                   m = partial.sd2.m.tt, 
                                   sd2 = partial.sd2.sd2.tt)
    partial.sd2.tt <- rbind( partial.sd2.tt,  partial.sd2.tt.i)
    
    
    partial.gamma.m.t1t <- exp(-1/gamma)*tt.record.i$m/gamma^2 + exp(-1/gamma)*partial.gamma.m.tt
    partial.gamma.sd2.t1t <- -2*exp(-2/gamma)/gamma^2+ exp(-2/gamma)*partial.gamma.sd2.tt + 2/gamma^2*exp(-2/gamma)*tt.record.i$sigma^2
    
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
    param.partial.m <- (y[i+1] - t1t.record.i$m)/t1t.record.i$sigma^2 
    param.partial.sd2 <- 0.5/t1t.record.i$sigma^2 - 0.5*(y[i+1]-t1t.record.i$m)^2/t1t.record.i$sigma^4
    
    logp.sd2 <-  param.partial.m*partial.sd2.i$m -param.partial.sd2*partial.sd2.i$sd2
    logp.gamma <-  param.partial.m*partial.gamma.i$m -param.partial.sd2*partial.gamma.i$sd2
    
    partial.logp.i <- data.frame(k = i, 
                                 partial.gamma = logp.gamma, 
                                 partial.sd2 = logp.sd2)
    partial.logp <- rbind(partial.logp, partial.logp.i)
  }  
  
  # compute sum of each log-likelihood term
  partial.gamma.logp <- sum(partial.logp$partial.gamma)
  partial.sd2.logp <- sum(partial.logp$partial.sd2) - 0.5/(1+sigma.w^2)
  
  
  return(list(gamma = partial.gamma.logp, sd2 = partial.sd2.logp))
}

# compute log-likelihood
compute_logp <- function(y, gamma, sigma.w, m0, sigma0){
  # compute the KF inference
  KF.result <- gp.KF(y, gamma, sigma.w, m0, sigma0)
  t1t.record <- KF.result$t1t
  
  # compute log-likelihood
  n <- length(y)
  y.mean <- t1t.record$m[-n]
  y.sd <- t1t.record$sigma[-n]+sigma.w^2
  value <- sum(log(dnorm(y[2:n], mean=y.mean, sd=y.sd))) + 
    log(dnorm(y[1], mean = 0, sd = 1+sigma.w^2))
  
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

