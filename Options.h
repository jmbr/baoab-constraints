#ifndef OPTIONS_H
#define OPTIONS_H

#include <cstdlib>

struct Options {
  double temperature;
  double friction;
  double min_dt, max_dt;
  unsigned num_experiments;
  double total_time;
  unsigned long random_seed;
  unsigned nbins;
  bool plot;

  Options()
      : temperature(-1),
        friction(1e4),
        min_dt(-1),
        max_dt(-1),
        num_experiments(0),
        total_time(-1),
        random_seed(0UL),
        nbins(50),
        plot(false) {}

  Options(double temperature_,
          double friction_,
          double min_dt_,
          double max_dt_,
          unsigned num_experiments_,
          double total_time_,
          unsigned long random_seed_,
          unsigned nbins_,
          bool plot_ = false)
      : temperature(temperature_),
        friction(friction_),
        min_dt(min_dt_),
        max_dt(max_dt_),
        num_experiments(num_experiments_),
        total_time(total_time_),
        random_seed(random_seed_),
        nbins(nbins_),
        plot(plot_) {}

  int parse(int argc, char* argv[]);
};

#endif
