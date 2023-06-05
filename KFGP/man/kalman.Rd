\name{kalman}
\alias{kalman}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inference of GP using Kalman Filter
}
\description{
The function returns a the predicted mean, predicted
variance, updated mean, updated variance, smoothing mean and smoothing
variance.
}
\usage{
kalman(y, phi, Sigmav, Sigmaw, m0, Sigma0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
a vector of state-space time series
}
  \item{phi}{
parameter, phi = exp(-1/gamma)
}
  \item{Sigmav}{
parameter, sigma_v^2 = 1-exp(-2/gamma)
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
