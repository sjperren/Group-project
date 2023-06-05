\name{optim_parm}
\alias{optim_parm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hyper-parameter selection
}
\description{
By applying the empirical Bayesian method, we optimised the hyperparameters (gamma, sigmaW) by finding their maximizers of log-likelihood.
}
\usage{
optim_parm(y, m0 = 0, Sigma0 = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
a vector of state-space time series
}
  \item{m0}{
parameter, initialised mean
}
  \item{Sigma0}{
parameter, initialised sd
}
}
