\name{kalmanSimulate}
\alias{kalmanSimulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a time series following a state-space model
}
\description{
simulate the state space model of the form
X_t = \phi X_{t-1} + V, Y_t = X_t + W
where
W\sim \mathcal{N}(0,\sigma_w^2)$ and $V\sim \mathcal{N}(0,\sigma_v^2).
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
