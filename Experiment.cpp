#include <cstdlib>
#include <cmath>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "Experiment.h"


Experiment::Experiment()
    : total_steps(0), plot(false) {}

Experiment::Experiment(double friction,
                       double temperature,
                       double dt_,
                       double total_time,
                       unsigned long random_seed,
                       bool plot_)
    : baoab(friction, temperature, dt_, random_seed),
      total_steps(static_cast<unsigned long long>(ceil(total_time / dt_))),
      plot(plot_),
      files_are_open(false) {}

Experiment::~Experiment() {
  if (files_are_open)
    closeFiles();
}

Experiment& Experiment::operator=(const Experiment& other) {
  if (this != &other) {
    baoab = other.baoab;
    average = other.average;
    total_steps = other.total_steps;
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
      << "total steps = " << total_steps << ", "
      << std::endl;

  if (plot) {
    plt1.open();
    plt2.open();
  }

  files_are_open = true;
}

void Experiment::closeFiles() {
  log.close();
  results.close();

  if (plot) {
    plt1.close();
    plt2.close();
  }

  files_are_open = false;
}

void Experiment::simulate() {
  for (size_t step = 1; step <= total_steps; step++) {
    baoab.advance();

    average.update(baoab.end_to_end_distance());

    if (step % static_cast<size_t>(1e4) == 0 || step == total_steps) {
      // if (plot) {
      //   baoab.plot(plt2);
      // }

      const double t = static_cast<double>(step) * baoab.dt;

      log << t << " "
          << trans(baoab.q) << " "
          << trans(baoab.p) << "\n\n"
          << std::endl;

      results << t << " " << average() << std::endl;
    }
  }
}
