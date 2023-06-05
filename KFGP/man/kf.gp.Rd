\name{kf.gp}
\alias{kf.gp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Smoothing and Srediction under Kalman Filter
}
\description{
This function embeds the Kalman filter algorithm into fitting Gaussian process of state-space time series
}
\usage{
kf.gp(y, gamma, Sigmaw, m0 = 0, Sigma0 = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
a vector of state-space time series
}
  \item{gamma}{
parameter, phi = exp(-1/gamma) and sigma_v^2 = 1-exp(-2/gamma)
}
  \item{Sigmaw}{
parameter, sd of noise term
}
  \item{m0}{
parameter, initialised mean
}
  \item{Sigma0}{
parameter, initialised sd
}
}
