#include <cstdlib>
#include <cmath>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "BAOAB.h"
#include "Average.h"
#include "Plotter.h"
#include "Experiment.h"


Experiment::Experiment()
    : equilibration_steps(0UL),
      production_steps(0ULL),
      plot(false) {}

Experiment::Experiment(double friction,
                       double temperature,
                       double dt_,
                       double equilibration_time,
                       double production_time,
                       unsigned long random_seed,
                       bool plot_)
    : baoab(friction, temperature, dt_, random_seed),
      equilibration_steps(long(ceil(equilibration_time / dt_))),
      production_steps(long(ceil(production_time / dt_))),
      plot(plot_),
      files_are_open(false) {}

Experiment::~Experiment() {
  if (files_are_open)
    closeFiles();
}

Experiment& Experiment::operator=(const Experiment& other) {
  if (this != &other) {
    baoab = other.baoab;
    end_to_end = other.end_to_end;
    potential = other.potential;
    equilibration_steps = other.equilibration_steps;
    production_steps = other.production_steps;
    plot = other.plot;
  }

  return *this;
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

void Experiment::simulate() {
  // Equilibrate
  for (size_t step = 1; step <= equilibration_steps; step++) {
    if (step % size_t(1e5) == 0) {
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

    // Collect ensemble averages.
    end_to_end.update(baoab.end_to_end_distance());
    potential.update(baoab.potential());

    if (step % size_t(1e4) == 0 || step == 1 || step == production_steps) {
      if (plot)
        baoab.plot(plt);

      const double t = double(step) * baoab.dt;

      log << t << "\n"
          << trans(baoab.q) << trans(baoab.p)
          << std::endl;

      results << t << " "
              << end_to_end.value() << " "
              << potential.value()
              << std::endl;
    }
  }
}
