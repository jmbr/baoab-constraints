#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <cstdlib>

#include <fstream>

class BAOAB;
class Average;
class Plotter;

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

 public:
  const unsigned long random_seed;
  BAOAB baoab;
  Average end_to_end;
  Average potential;
  unsigned long equilibration_steps;
  unsigned long long production_steps;
  const bool plot;
  Plotter plt;

 private:
  bool files_are_open;
  std::ofstream log;
  std::ofstream results;
};

#endif
