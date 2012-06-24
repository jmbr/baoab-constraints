#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <cstdlib>

#include <fstream>

#include "BAOAB.h"
#include "Histogram.h"
#include "Plotter.h"

class Experiment {
 public:
  Experiment();

  Experiment(double K, double friction, double temperature,
             double dt, double time,
             unsigned random_seed, unsigned nbins,
             bool plot = false);

  Experiment(const Experiment& e);

  ~Experiment();

  Experiment& operator=(const Experiment& other) {
    if (this != &other) {
      baoab = other.baoab;
      histogram = other.histogram;
      total_steps = other.total_steps;
      plot = other.plot;
    }

    return *this;
  }

  void openFiles();
  void closeFiles();

  void simulate();

 private:
  void compute_step();

 public:
  BAOAB_with_RATTLE baoab;
  // BAOAB_with_DoPri baoab;
  // BAOAB_with_Rotation baoab;
  Histogram histogram;
  unsigned long long total_steps;
  bool plot;
  Plotter plt1, plt2;

 private:
  bool files_are_open;
  std::ofstream log;
  std::ofstream results;
};

#endif
