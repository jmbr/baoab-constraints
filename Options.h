#ifndef OPTIONS_H
#define OPTIONS_H

#include <cstdlib>

struct Options {
  double temperature;
  double friction;
  double min_dt, max_dt;
  unsigned num_experiments;
  double equilibration_time;
  double production_time;
  unsigned long random_seed;
  bool plot;

  Options()
      : temperature(-1),
        friction(1e4),
        min_dt(-1),
        max_dt(-1),
        num_experiments(0),
        equilibration_time(1e7),
        production_time(-1),
        random_seed(0UL),
        plot(false) {}

  Options(double temperature_,
          double friction_,
          double min_dt_,
          double max_dt_,
          unsigned num_experiments_,
          double equilibration_time_,
          double production_time_,
          unsigned long random_seed_,
          bool plot_ = false)
      : temperature(temperature_),
        friction(friction_),
        min_dt(min_dt_),
        max_dt(max_dt_),
        num_experiments(num_experiments_),
        equilibration_time(equilibration_time_),
        production_time(production_time_),
        random_seed(random_seed_),
        plot(plot_) {}

  int parse(int argc, char* argv[]);
};

#endif
