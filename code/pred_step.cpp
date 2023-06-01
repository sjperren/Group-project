#include <Rcpp.h>
using namespace Rcpp;

// C++ function to perform the prediction step
List pred_step(double m_tt, double sigma_tt, double phi, double sigma_v) {
  double m_t1t = phi * m_tt;
  double sigma_t1t = sqrt(sigma_v * sigma_v + phi * phi * sigma_tt * sigma_tt);
  
  return List::create(Named("m") = m_t1t, Named("sigma") = sigma_t1t);
}

// Rcpp module to export the C++ function
RCPP_MODULE(MyModule) {
  function("pred_step", &pred_step, "Perform the prediction step");
}

// Define the entry point to the DLL
extern "C" void R_init_MyPackage(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
