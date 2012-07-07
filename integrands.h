#ifndef INTEGRANDS_H
#define INTEGRANDS_H

#ifdef __cplusplus
extern "C" {
#endif

  double end_to_end_gsl(double x[], size_t dim, void* param);
  void end_to_end_cub(unsigned ndim, const double x[],
                      void* fdata, unsigned fdim, double* fval);

  double partition_function_gsl(double x[], size_t dim, void* param);
  void partition_function_cub(unsigned ndim, const double x[],
                              void* fdata, unsigned fdim, double* fval);

#ifdef __cplusplus
}
#endif

#endif
