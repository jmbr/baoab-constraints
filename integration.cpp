#include <cmath>
#include <cassert>
#include <cstdlib>

#include <iomanip>
#include <iostream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte_vegas.h>

extern "C" {
#include "integrands.h"
}

const double kT = 0.5;
const size_t dim = 6;
const size_t num_evals = size_t(1e7);

int main() {
  gsl_rng_env_setup();
  
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
  
  gsl_monte_vegas_state* state = gsl_monte_vegas_alloc(dim);
  gsl_monte_vegas_init(state);

  gsl_monte_function f1 = { &partition_function_integrand, dim, (void*) &kT };
  gsl_monte_function f2 = { &end_to_end_integrand, dim, (void*) &kT };

  double xl[] = { -M_PI, -M_PI, -M_PI, -M_PI, -M_PI, -M_PI };
  double xu[] = { +M_PI, +M_PI, +M_PI, +M_PI, +M_PI, +M_PI };

  double partition_function, ensemble_average, abserr;
  
  gsl_monte_vegas_integrate(&f1, xl, xu, dim, num_evals,
                            rng, state, &partition_function, &abserr);

  std::cout << std::setprecision(10)
            << "Partition function: " << partition_function
            << " (abs. err.: " << abserr << ")"
            << std::endl;

  gsl_monte_vegas_integrate(&f2, xl, xu, dim, num_evals,
                            rng, state, &ensemble_average, &abserr);

  ensemble_average /= partition_function;
  
  std::cout << std::setprecision(10)
            << "Ensemble average: " << ensemble_average
            << " (abs. err.: " << abserr << ")"
            << std::endl;

  gsl_monte_vegas_free(state);
  gsl_rng_free(rng);

  return 0;
}
