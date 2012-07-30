#ifndef INTEGRANDS_H
#define INTEGRANDS_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

  double end_to_end_integrand(double angles[], size_t dim, double kT);

  double partition_function_integrand(double angles[], size_t dim, double kT);

#ifdef __cplusplus
}
#endif

#endif
