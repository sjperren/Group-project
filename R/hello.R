# Exponential kernel
kernel <- function(a, b, gamma) {
  exp(-abs(a - b) / gamma)
}

# Covariance matrix
cov_matrix <- function(time_points, gamma) {
  n <- length(time_points)
  K <- matrix(nrow=n, ncol=n)

  for (i in 1:n) {
    for (j in 1:n) {
      K[i, j] <- kernel(time_points[i], time_points[j], gamma)
    }
  }

  return(K)
}

# Gaussian Process Regression
gp_regression <- function(time_points, observations, new_time_point, gamma) {
  K <- cov_matrix(time_points, gamma)
  k <- sapply(time_points, function(tp) kernel(tp, new_time_point, gamma))

  K_inv <- solve(K)

  mean <- sum(k %*% K_inv %*% observations)
  variance <- kernel(new_time_point, new_time_point, gamma) - sum(k %*% K_inv %*% k)

  return(list(mean = mean, variance = variance))
}

# Time points
t <- c(1, 2, 3, 4, 5)
# Observed values
y <- c(1, 2, 3, 2, 1)
# Gamma parameter
gamma <- 1

# New time point to predict
new_t <- 6

# Prediction
prediction <- gp_regression(t, y, new_t, gamma)
print(prediction)

