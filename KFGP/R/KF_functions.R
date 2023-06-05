
# Simulation of state-space model
kalmanSimulate <- function(T=100,m0=0,Sigma0=1,phi=1,Sigmav=0.02,Sigmaw=0.2){
  x <- rep(NA,T)
  y<- rep(NA,T)

  #create the hidden states X_t
  x[1]<-rnorm(1,mean=m0,sd=sqrt(Sigma0))
  for (i in 2:T){
    x[i] <- phi*x[i-1] + rnorm(1,mean=0,sd=sqrt(Sigmav))
  }

  #create the observed states Y_t
  for (i in 1:T){
    y[i] <- x[i] + rnorm(1,mean=0,sd=sqrt(Sigmaw))
  }
  return (list("x"=x,"y"=y))
}

# ----------------------------------------------------------------------------------------------------------------------
# Kalman filter
kalman <- function(y, phi, Sigmav, Sigmaw, m0, Sigma0){

  T <- length(y)

  #initialization
  mu.p <- rep(NA,T)
  Sigma.p <- rep(NA,T)
  mu.f <- rep(NA,T)
  Sigma.f <- rep(NA,T)
  mu.s <- rep(NA,T)
  Sigma.s <- rep(NA,T)


  #forward recursion time1
  mu.p[1] <- m0
  Sigma.p[1] <- Sigma0
  mu.f[1] <- m0 + (y[1]-m0)*(Sigma0/(Sigma0+Sigmaw))
  Sigma.f[1] <- Sigma0-(Sigma0^2/(Sigma0+Sigmaw))

  #forward recursion time 2:T
  for (t in 2:T){

    #prediction
    mu.p[t] <- phi*mu.f[t-1]
    Sigma.p[t] <- phi^2 * Sigma.f[t-1] + Sigmaw

    #update
    deno <- Sigmaw + Sigma.p[t]
    mu.f[t] <- Sigmaw*mu.p[t]/deno + Sigma.p[t]*y[t]/deno
    Sigma.f[t] <- Sigmaw*Sigma.p[t]/deno
  }

  #backword recursion
  mu.s[T] <- mu.f[T]
  Sigma.s[T] <- Sigma.f[T]
  for (t in (T-1):1){
    J <- phi*Sigma.f[t]/Sigma.p[t+1]
    mu.s[t] <- mu.f[t] + J*(mu.s[t+1]-mu.p[t+1])
    Sigma.s[t] <- Sigma.f[t] + J^2*(Sigma.s[t+1]-Sigma.p[t+1])
  }
  return (list(mu.f=mu.f, Sigma.f=Sigma.f,
               mu.p=mu.p, Sigma.p=Sigma.p,
               mu.s=mu.s, Sigma.s=Sigma.s))
}

# ----------------------------------------------------------------------------------------------------------------------
# smoothing and prediction under KF
kf.gp <- function(y,gamma,Sigmaw,m0=0,Sigma0=1){
  T <- length(y)
  #update Sigmav and phi
  phi <- exp(-1/gamma)
  Sigmav <- 1-exp(-2/gamma)
  result <- kalman(y, phi=phi,
                   Sigmav=Sigmav, Sigmaw=Sigmaw,
                   m0=m0, Sigma0=Sigma0)

  return (list(mu.p=result$mu.p, Sigma.p=result$Sigma.p,
               mu.f=result$mu.f, Sigma.f=result$Sigma.f,
               mu.s=result$mu.s, Sigma.s=result$Sigma.s))

}
# ----------------------------------------------------------------------------------------------------------------------
# compute log-likelihood function
kf.loglikelihood1 <- function(y,mu.p,Sigma.p,Sigmaw,m0=0,Sigma0=1){
  T <- length(y)
  likelihood <- rep(NA,T)

  #at time 1
  likelihood[1] <- log(dnorm(y[1], mean = m0,
                             sd = sqrt(Sigma0 + Sigmaw)))

  #time 2:T
  for (t in 2:T){
    likelihood[t] <- log(dnorm(y[t], mean = mu.p[t],
                               sd = sqrt(Sigmaw+Sigma.p[t])))
  }
  return (sum(likelihood))
}

# ----------------------------------------------------------------------------------------------------------------------
# Optimize over Loglikelihood
kf.loglikelihood <- function(y,gamma,Sigmaw,m0=0,Sigma0=1){
  o <- kf.gp(y=y,gamma=gamma,Sigmaw=Sigmaw,m0=m0,Sigma0=Sigma0)
  mu.p <- o$mu.p
  Sigma.p <- o$Sigma.p
  result <- kf.loglikelihood1(y=y,mu.p=mu.p,Sigma.p=Sigma.p,Sigmaw=Sigmaw,m0=m0,Sigma0=Sigma0)
  return (result)
}

optim_parm <- function(y, m0=0, Sigma0=1){
  opt_param <- optim(par = c(5,0.5),
                     fn = function(parm) -1*kf.loglikelihood(y,parm[1],  parm[2], m0, Sigma0))
  return(list(gamma = opt_param$par[1],
              Sigmaw=opt_param$par[2]))
}

