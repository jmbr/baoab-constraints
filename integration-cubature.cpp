#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <iomanip>
#include <iostream>

#include <cubature.h>

#include "integrands.h"

using namespace std;

int main(int argc, char* argv[]) {
  if (argc < 3) {
    cerr << "Usage: " << argv[0] << " TEMPERATURE EVALUATIONS"
              << endl;
    return EXIT_FAILURE;
  }

  const double kT = atof(argv[1]);
  const unsigned num_evals = unsigned(atof(argv[2]));

  cout << setprecision(10)
            << "Estimating ensemble average at kT = " << kT
            << " by evaluating the integrand 10^"
            << log10(num_evals) << " times."
            << endl;

  const double xl[] = { -M_PI, -M_PI, -M_PI, -M_PI, -M_PI };
  const double xu[] = { +M_PI, +M_PI, +M_PI, +M_PI, +M_PI };
  const unsigned dim = sizeof(xl) / sizeof(xl[0]);

  double partition_function, partition_function_err;
  {
    int status;
    status = adapt_integrate(1, partition_function_cub, (void*) &kT,
                             dim, xl, xu, num_evals, 1e-8, 1e-10,
                             &partition_function, &partition_function_err);
    if (status != 0) {
      perror("adapt_integrate");
      return EXIT_FAILURE;
    }
  
    cout << "Partition function: " << partition_function
         << " (estimated error: " << partition_function_err << ")"
         << endl;
  }

  double ensemble_average, ensemble_average_err;
  {
    int status;
    status = adapt_integrate(1, end_to_end_cub, (void*) &kT,
                             dim, xl, xu, num_evals, 1e-8, 1e-10,
                             &ensemble_average, &ensemble_average_err);
    if (status != 0) {
      perror("adapt_integrate");
      return EXIT_FAILURE;
    }

    ensemble_average /= partition_function;
    cout << "Ensemble average: " << ensemble_average
         << " (estimated error: " << ensemble_average_err << ")"
         << endl;
  }
  
  return EXIT_SUCCESS;
}
