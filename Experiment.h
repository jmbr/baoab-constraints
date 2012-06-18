#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <cstdlib>

#include <fstream>

#include "BAOAB.h"
#include "Histogram.h"

class Experiment {
 public:
  Experiment(double K, double friction, double temperature,
             double dt, double time,
             unsigned random_seed, size_t nbins,
             bool plot = false);

  ~Experiment();
  
  void simulate();
  
 private:
  void compute_step();

 public:
  BAOAB_with_RATTLE baoab;
  // BAOAB_with_DoPri baoab;
  // BAOAB_with_Rotation baoab;
  Histogram histogram;
  double dt;
  unsigned long long total_steps;
  bool plot;
  
 private:
  std::ofstream log;
  std::ofstream results;
};

#endif
