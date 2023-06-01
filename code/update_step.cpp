#include <Rcpp.h>
using namespace Rcpp;

// C++ function to perform the update step
List update_step(double m_t1t, double sigma_t1t, double phi, double sigma_w, double new_y) {
  double denominator = sigma_w * sigma_w + sigma_t1t * sigma_t1t;
  
  double m_t1t1 = (sigma_w * sigma_w * m_t1t + sigma_t1t * new_y) / denominator;
  double sigma_t1t1 = sqrt(sigma_w * sigma_w * sigma_t1t * sigma_t1t / denominator);
  
  return List::create(Named("m") = m_t1t1, Named("sigma") = sigma_t1t1);
}

// Rcpp module to export the C++ function
RCPP_MODULE(MyModule) {
  function("update_step", &update_step, "Perform the update step");
}

// Define the entry point to the DLL
extern "C" void R_init_MyPackage(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
