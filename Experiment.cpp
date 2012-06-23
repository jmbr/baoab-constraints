#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>
#include <sstream>

#include "Experiment.h"


Experiment::Experiment()
    : total_steps(0), plot(false) {}

Experiment::Experiment(double K, double friction, double temperature,
                       double dt_, double total_time,
                       unsigned random_seed, unsigned nbins, bool plot_)
    : baoab(K, friction, temperature, dt_, random_seed),
      histogram(nbins, temperature, K),
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
    histogram = other.histogram;
    total_steps = other.total_steps;
    plot = other.plot;
  }

  return *this;
}

void Experiment::openFiles() {
  const double dt = baoab.dt;

  std::ostringstream log_name;
  log_name << "log-dt-" << dt << ".dat";
  log.open(log_name.str());
  if (!log.is_open()) {
    std::cerr << "Unable to open log file for dt = " << dt
              << std::endl;
    abort();                            // XXX Exception.
  }

  log << std::setprecision(14);

  std::ostringstream results_name;
  results_name << "result-dt-" << dt << ".dat";
  results.open(results_name.str());
  if (!results.is_open()) {
    std::cerr << "Unable to open results file for dt = "
              << dt << std::endl;
    abort();                            // XXX Exception.
  }

  results << std::setprecision(14);

  log << "K = " << baoab.K << ", "
      << "temperature = " << baoab.temperature << ", "
      << "friction = " << baoab.friction << ", "
      << "time step length = " << dt << ", "
      << "total steps = " << total_steps << ", "
      << "number of bins in histogram = " << histogram.size()
      << std::endl;

  plt1.open();
  plt2.open();

  files_are_open = true;
}

void Experiment::closeFiles() {
  log.close();
  results.close();

  plt1.close();
  plt2.close();

  files_are_open = false;
}

void Experiment::compute_step() {
  baoab();
  ++histogram[baoab.angle()];
}

void Experiment::simulate() {
  for (size_t step = 1; step <= total_steps; step++) {
    compute_step();

    if (step % static_cast<size_t>(1e6) == 0 || step == total_steps) {
      if (plot) {
        histogram.plot();
        // baoab.plot();
      }

      const double t = static_cast<double>(step) * baoab.dt;

      log << t << " " << baoab.q << " " << baoab.p
          << "\n\n"
          << histogram << std::endl;

      results << t << " " << histogram.error() << std::endl;
    }
  }
}
