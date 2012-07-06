#ifndef INTEGRANDS_H
#define INTEGRANDS_H

#include <stdlib.h>

double partition_function_integrand(double angles[], size_t dim, void* params);

double end_to_end_integrand(double angles[], size_t dim, void* params);

#endif
