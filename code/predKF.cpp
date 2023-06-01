#include <Rcpp.h>
using namespace Rcpp;

// C++ function to perform prediction in Kalman Filter
List predKF(double m_nn, double sd_nn, double gamma) {
  // Compute the parameters
  double phi = exp(-1 / gamma);
  double sigma_v = sqrt(1 - exp(-2 / gamma));
  
  // Predicted m
  double pred_m = phi * m_nn;
  
  // Predicted sd^2
  double pred_sd2 = pow(sd_nn, 2) * pow(phi, 2) + pow(sigma_v, 2);
  
  // Output
  return List::create(Named("m") = pred_m, Named("sigma") = sqrt(pred_sd2));
}

// Rcpp module to export the C++ function
RCPP_MODULE(MyModule) {
  function("predKF", &predKF, "Perform prediction in Kalman Filter");
}

// Define the entry point to the DLL
extern "C" void R_init_MyPackage(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
