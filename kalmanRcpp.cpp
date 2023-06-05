#include <Rcpp.h>
using namespace Rcpp;

// Rcpp version of the function
// The function takes NumericVector y, double phi, double Sigmav, double Sigmaw, double m0, and double Sigma0 as input
// It returns a list containing the calculated values
List kalmanRcpp(NumericVector y, double phi, double Sigmav, double Sigmaw, double m0, double Sigma0) {
  int T = y.size();
  
  // Initialization
  NumericVector mu_p(T);
  NumericVector Sigma_p(T);
  NumericVector mu_f(T);
  NumericVector Sigma_f(T);
  NumericVector mu_s(T);
  NumericVector Sigma_s(T);
  
  // Forward recursion time 1
  mu_p[0] = m0;
  Sigma_p[0] = Sigma0;
  mu_f[0] = m0 + (y[0] - m0) * (Sigma0 / (Sigma0 + Sigmaw));
  Sigma_f[0] = Sigma0 - pow(Sigma0, 2) / (Sigma0 + Sigmaw);
  
  // Forward recursion time 2:T
  for (int t = 1; t < T; t++) {
    // Prediction
    mu_p[t] = phi * mu_f[t - 1];
    Sigma_p[t] = pow(phi, 2) * Sigma_f[t - 1] + Sigmaw;
    
    // Update
    double deno = Sigmaw + Sigma_p[t];
    mu_f[t] = (Sigmaw * mu_p[t]) / deno + (Sigma_p[t] * y[t]) / deno;
    Sigma_f[t] = (Sigmaw * Sigma_p[t]) / deno;
  }
  
  // Backward recursion
  mu_s[T - 1] = mu_f[T - 1];
  Sigma_s[T - 1] = Sigma_f[T - 1];
  for (int t = T - 2; t >= 0; t--) {
    double J = phi * Sigma_f[t] / Sigma_p[t + 1];
    mu_s[t] = mu_f[t] + J * (mu_s[t + 1] - mu_p[t + 1]);
    Sigma_s[t] = Sigma_f[t] + pow(J, 2) * (Sigma_s[t + 1] - Sigma_p[t + 1]);
  }
  
  // Create and return the list
  List result = List::create(
    Named("mu.f") = mu_f,
    Named("Sigma.f") = Sigma_f,
    Named("mu.p") = mu_p,
    Named("Sigma.p") = Sigma_p,
    Named("mu.s") = mu_s,
    Named("Sigma.s") = Sigma_s
  );
  
  return result;
}

// Rcpp module to expose the function to R
RCPP_MODULE(kalmanModule) {
  function("kalmanRcpp", &kalmanRcpp);
}
