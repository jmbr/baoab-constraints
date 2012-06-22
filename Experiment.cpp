#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>
#include <sstream>

#include "Experiment.h"

Experiment::Experiment(double K, double friction, double temperature,
                       double dt_, double total_time,
                       unsigned random_seed, size_t nbins, bool plot_)
    : baoab(K, friction, temperature, dt_, random_seed),
      histogram(nbins, temperature, K),
      dt(dt_),
      plot(plot_) {
  total_steps = static_cast<unsigned long long>(std::ceil(total_time / dt));

  std::ostringstream log_name;
  log_name << "log-dt-" << dt << ".dat";
  log.open(log_name.str());
  if (!log.is_open()) {
    std::cerr << "Unable to open log file for dt = " << dt << std::endl;
    abort();                            // XXX Exception.
  }

  log << std::setprecision(14);

  std::ostringstream results_name;
  results_name << "result-dt-" << dt << ".dat";
  results.open(results_name.str());
  if (!results.is_open()) {
    std::cerr << "Unable to open results file for dt = " << dt << std::endl;
    abort();                            // XXX Exception.
  }

  results << std::setprecision(14);

  log << "K = " << K << ", "
      << "temperature = " << temperature << ", "
      << "friction = " << friction << ", "
      << "time step length = " << dt << ", "
      << "total time = " << total_time << ", "
      << "random seed = " << random_seed << ", "
      << "number of bins in histogram = " << nbins
      << std::endl;
}

Experiment::~Experiment() {
  log.close();
  results.close();
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

      const double t = static_cast<double>(step) * dt;

      log << t << " " << baoab.q << " " << baoab.p
          << "\n\n"
          << histogram << std::endl;

      results << t << " " << histogram.error() << std::endl;
    }
  }
}
