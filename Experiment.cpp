#include <cstdlib>
#include <cmath>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "BAOAB.h"
#include "Average.h"
#include "Plotter.h"
#include "Histogram.h"
#include "Experiment.h"

Experiment::Experiment(double friction,
                       double temperature,
                       double dt_,
                       double equilibration_time,
                       double production_time,
                       unsigned long random_seed_,
                       bool plot_)
    : random_seed(random_seed_),
      baoab(friction, temperature, dt_, random_seed),
      equilibration_steps(long(ceil(equilibration_time / dt_))),
      production_steps(long(ceil(production_time / dt_))),
      plot(plot_),
      histogram(40, temperature, K),
      files_are_open(false),
      update_interval(1e5) {
}

Experiment::~Experiment() {
  if (files_are_open)
    closeFiles();
}

void Experiment::openFiles() {
  const double dt = baoab.dt;

  std::ostringstream log_name;
  log_name << "log-dt-" << dt << ".dat";
  log.open(log_name.str().c_str());
  if (!log.is_open()) {
    std::cerr << "Unable to open log file for dt = " << dt
              << std::endl;
    abort();                            // XXX Exception.
  }

  log << std::setprecision(14);

  std::ostringstream results_name;
  results_name << "result-dt-" << dt << ".dat";
  results.open(results_name.str().c_str());
  if (!results.is_open()) {
    std::cerr << "Unable to open results file for dt = "
              << dt << std::endl;
    abort();                            // XXX Exception.
  }

  results << std::setprecision(14);

  log << "temperature = " << baoab.temperature << ", "
      << "friction = " << baoab.friction << ", "
      << "time step length = " << dt << ", "
      << "equilibration steps = " << equilibration_steps << ", "
      << "production steps = " << production_steps << ", "
      << "random seed = " << random_seed
      << std::endl;

  if (plot)
    plt.open();

  files_are_open = true;
}

void Experiment::closeFiles() {
  log.close();
  results.close();

  if (plot)
    plt.close();

  files_are_open = false;
}

bool Experiment::should_update(size_t step) const {
  return step == 1
      || step % update_interval == 0
      || step == production_steps;
}

void Experiment::simulate() {
  // Equilibrate
  for (size_t step = 1; step <= equilibration_steps; step++) {
    if (step % update_interval == 0) {
      log << "Running equilibraton phase ("
          << floor(double(step) / double(equilibration_steps) * 100.0)
          << "% completed)."
          << std::endl;
    }

    baoab.advance();
  }

  // Do production simulation
  log << "Running production simulation." << std::endl;
  for (size_t step = 1; step <= production_steps; step++) {
    baoab.advance();
    ++histogram[baoab.angle()];

    // Collect ensemble averages.
    potential.update(baoab.potential());
    // end_to_end.update(baoab.end_to_end_distance());

    if (should_update(step)) {
      baoab.center(); // XXX Make sure this is not messing up the averages.

      if (plot)
        histogram.plot(plt);
        // baoab.plot(plt);

      const double t = double(step) * baoab.dt;

      log << t << "\n"
          << trans(baoab.q) << trans(baoab.p)
          << std::endl;

      results << t << " "
              // << end_to_end.value() << " "
              << potential.value() << " "
              << histogram.error()
              << std::endl;
    }
  }
}
