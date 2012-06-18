#ifndef OPTIONS_H
#define OPTIONS_H

#include <cstdlib>

struct Options {
  double K;
  double temperature;
  double friction;
  double min_dt, max_dt;
  size_t num_experiments;
  double total_time;
  unsigned random_seed;
  size_t nbins;
  bool plot;

  Options()
      : K(0.5),
        temperature(0.1),
        friction(1e4),
        min_dt(0.35),
        max_dt(0.43),
        num_experiments(21),
        total_time(1e10),
        random_seed(0),
        nbins(50),
        plot(false) {}
  
  Options(double K_,
          double temperature_,
          double friction_,
          double min_dt_,
          double max_dt_,
          size_t num_experiments_,
          double total_time_,
          unsigned random_seed_,
          size_t nbins_,
          bool plot_ = false)
      : K(K_),
        temperature(temperature_),
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
