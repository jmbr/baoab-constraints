#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <iostream>

#include <nag.h>
#include <nag_stdlib.h>
#include <nagd01.h>

#include "integrands.h"

extern "C" {
  static double NAG_CALL partfun(Integer n, double z[], Nag_User *comm);

  static double NAG_CALL ensavg(Integer n, double z[], Nag_User *comm);
}

using namespace std;

int main(int argc, char *argv[]) {
  if (argc < 3) {
    cerr << "Usage: " << argv[0] << " TEMPERATURE EVALUATIONS"
         << endl;
    return EXIT_FAILURE;
  }

  double kT = atof(argv[1]);
  Integer maxpts = Integer(atof(argv[2]));

  double a[] = { -M_PI, -M_PI, -M_PI, -M_PI, -M_PI, -M_PI };
  double b[] = { +M_PI, +M_PI, +M_PI, +M_PI, +M_PI, +M_PI };
  const Integer ndim = sizeof(a) / sizeof(a[0]);

  Nag_User comm;
  comm.p = &kT;

  NagError fail;
  INIT_FAIL(fail);

  Nag_MCMethod method = Nag_ManyIterations;
  Nag_Start    cont = Nag_Cold;

  double partfun_finval, ensavg_finval;
  double partfun_acc, ensavg_acc;
  Integer partfun_minpts = 0, ensavg_minpts = 0;
  const double eps = 1e-12;

  fprintf(stderr, "Dimensionality = %ld\nRequested accuracy =%12.2e\n",
          ndim, eps);

  {
    double* comm_arr = 0;

    nag_multid_quad_monte_carlo_1(ndim, partfun, method, cont, a, b, &partfun_minpts, maxpts, eps, &partfun_finval, &partfun_acc, &comm_arr, &comm, &fail);

    if (fail.code != NE_NOERROR
        && fail.code != NE_QUAD_MAX_INTEGRAND_EVAL) {
      fprintf(stderr, "Error: nag_multid_quad_adapt_1: %s\n",
              fail.message);
      return EXIT_FAILURE;      
    }

    fprintf(stderr, "Estimated value    =%12.2e\n", partfun_finval);
    fprintf(stderr, "Estimated accuracy =%12.2e\n", partfun_acc);
    fprintf(stderr, "Number of evaluations = %5"NAG_IFMT"\n", partfun_minpts);

    if (comm_arr) NAG_FREE(comm_arr);
  }

  {
    double* comm_arr = 0;

    nag_multid_quad_monte_carlo_1(ndim, ensavg, method, cont, a, b, &ensavg_minpts, maxpts, eps, &ensavg_finval, &ensavg_acc, &comm_arr, &comm, &fail);

    if (fail.code != NE_NOERROR
        && fail.code != NE_QUAD_MAX_INTEGRAND_EVAL) {
      fprintf(stderr, "Error: nag_multid_quad_montecarlo_1: %s\n",
              fail.message);
      return EXIT_FAILURE;      
    }

    fprintf(stderr, "Estimated value    =%12.2e\n", ensavg_finval);
    fprintf(stderr, "Estimated accuracy =%12.2e\n", ensavg_acc);
    fprintf(stderr, "Number of evaluations = %5"NAG_IFMT"\n", ensavg_minpts);

    if (comm_arr) NAG_FREE(comm_arr);
  }

  printf("%ld %g %g %g\n",
         min(partfun_minpts, ensavg_minpts),
         ensavg_finval / partfun_finval,
         (ensavg_finval - ensavg_acc) / (partfun_finval + partfun_acc),
         (ensavg_finval + ensavg_acc) / (partfun_finval - partfun_acc));
  
  return EXIT_SUCCESS;
}

double NAG_CALL partfun(Integer n, double z[], Nag_User *comm) {
  const double kT = *((double*) comm->p);
  return partition_function_integrand(z, n, kT);
}

double NAG_CALL ensavg(Integer n, double z[], Nag_User *comm) {
  const double kT = *((double*) comm->p);
  return end_to_end_integrand(z, n, kT);
}
