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
             double time,
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
  Average average;
  unsigned long long total_steps;
  bool plot;
  Plotter plt;

 private:
  bool files_are_open;
  std::ofstream log;
  std::ofstream results;
};

#endif
