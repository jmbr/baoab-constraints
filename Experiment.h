#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <cstdlib>

#include <fstream>

class BAOAB;
class Average;
class Plotter;

class Experiment {
 public:
  Experiment();

  Experiment(double friction,
             double temperature,
             double dt,
             double equilibration_time,
             double production_time,
             unsigned long random_seed,
             bool plot = false);

  Experiment(const Experiment& e);

  ~Experiment();

  Experiment& operator=(const Experiment& other);

  void openFiles();
  void closeFiles();

  void simulate();

 private:
  void advance();

 public:
  BAOAB baoab;
  Average end_to_end;
  Average potential;
  unsigned long equilibration_steps;
  unsigned long long production_steps;
  bool plot;
  Plotter plt;
  unsigned long random_seed;

 private:
  bool files_are_open;
  std::ofstream log;
  std::ofstream results;
};

#endif
