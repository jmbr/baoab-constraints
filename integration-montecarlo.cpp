#include <cmath>
#include <cstdlib>

#include <iomanip>
#include <iostream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte_vegas.h>

#include "integrands.h"

using namespace std;

int main(int argc, char* argv[]) {
  if (argc < 3) {
    cerr << "Usage: " << argv[0] << " TEMPERATURE EVALUATIONS"
              << endl;
    return EXIT_FAILURE;
  }

  const double kT = atof(argv[1]);
  const size_t num_evals = size_t(atof(argv[2]));

  cout << setprecision(10)
            << "Estimating ensemble average at kT = " << kT
            << " by evaluating the integrand " << num_evals << " times."
            << endl;

  double xl[] = { -M_PI, -M_PI, -M_PI, -M_PI, -M_PI };
  double xu[] = { +M_PI, +M_PI, +M_PI, +M_PI, +M_PI };
  unsigned dim = sizeof(xl) / sizeof(xl[0]);
  
  gsl_rng_env_setup();
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
  
  gsl_monte_vegas_state* state = gsl_monte_vegas_alloc(dim);
  gsl_monte_vegas_init(state);

  gsl_monte_vegas_params params;
  gsl_monte_vegas_params_get(state, &params);
  params.iterations = 10;
  params.verbose = 0;
  gsl_monte_vegas_params_set(state, &params);
  
  double partition_function, partition_function_err;
  {
    gsl_monte_function f1 = { &partition_function_gsl, dim, (void*) &kT };
    gsl_monte_vegas_integrate(&f1, xl, xu, dim, num_evals, rng, state,
                              &partition_function, &partition_function_err);

    cout << "Partition function: " << partition_function
         << " (estimated error: " << partition_function_err << ")"
         << endl;
  }

  double ensemble_average, ensemble_average_err;
  {
    gsl_monte_function f2 = { &end_to_end_gsl, dim, (void*) &kT };
    gsl_monte_vegas_integrate(&f2, xl, xu, dim, num_evals, rng, state,
                              &ensemble_average, &ensemble_average_err);
    ensemble_average /= partition_function;

    cout << "Ensemble average: " << ensemble_average
         << " (estimated error: " << ensemble_average_err << ")"
         << endl;
  }

  gsl_monte_vegas_free(state);
  gsl_rng_free(rng);

  return 0;
}
