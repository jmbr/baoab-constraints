#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <cstdlib>

#include <fstream>

class BAOAB;
class Average;
class Plotter;
class Histogram;

class Experiment {
 public:
  Experiment(double friction,
             double temperature,
             double dt,
             double equilibration_time,
             double production_time,
             unsigned long random_seed,
             bool plot = false);

  ~Experiment();

  void openFiles();
  void closeFiles();

  void simulate();

 private:
  void advance();
  bool should_update(size_t step) const;

 public:
  const unsigned long random_seed;
  BAOAB baoab;
  Average end_to_end;
  Average potential;
  unsigned long equilibration_steps;
  unsigned long long production_steps;
  const bool plot;
  Plotter plt;
  Histogram histogram;

 private:
  bool files_are_open;
  std::ofstream log;
  std::ofstream results;
  const size_t update_interval;
};

#endif
