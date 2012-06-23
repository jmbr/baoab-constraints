#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <numeric>
#include <algorithm>
#include <functional>

#include "linspace.h"
#include "Histogram.h"
#include "Plotter.h"

using std::vector;

Histogram::Histogram(unsigned nbins_, double temperature_, double K_)
    : nbins(nbins_), histogram(nbins), reference(nbins),
      temperature(temperature_), K(K_) {
  compute();
}

static double f(double theta, void* params) {
  const double temperature = *((double *) params);
  const double K = *((double *) params + 1);

  double s, c;
  sincos(theta, &s, &c);
  const double D = sqrt(K) * c * c + s * s / sqrt(K);

  const double r = sqrt(c * c / K + s * s);
  const double d = 1.0, a = 2.0, r0 = 1.0;
  const double exp_factor = (1.0 - exp(-a * (r - r0)));
  const double U = d * exp_factor * exp_factor;

  // const double U = c * c * s / K;
  // const double U = c * c * s * s / K;
  // const double U = s;

  return sqrt(D) * exp(-U / temperature);
}

void Histogram::compute() {
  gsl_integration_workspace* w;

  w = gsl_integration_workspace_alloc(10000);

  double params[] = { temperature, K };

  gsl_function F;
  F.function = &f;
  F.params = &params[0];

  double err;
  double partition_function;
  vector<double> intervals = linspace<double>(-M_PI, M_PI, nbins+1);

  gsl_integration_qags(&F, -M_PI, M_PI, 0, 1e-12, 10000,
                       w, &partition_function, &err);

  for (unsigned n = 0; n < histogram.size(); n++) {
    const double a = intervals[n], b = intervals[n+1];

    gsl_integration_qags(&F, a, b, 0, 1e-12, 10000,
                         w, &reference[n], &err);

    reference[n] /= partition_function;
  }

  minimum = *std::min_element(reference.begin(), reference.end());
  maximum = *std::max_element(reference.begin(), reference.end());

#if 0
  std::cout << std::setprecision(16)
            << "Partition function = "
            << partition_function << "\n";
  for (unsigned k = 0; k < nbins; k++)
    std::cout << reference[k] << " ";
  std::cout << endl;
#endif

  gsl_integration_workspace_free(w);
}

unsigned long& Histogram::operator[](double theta) {
  const double a = -M_PI, b = M_PI;
  const double N = static_cast<double>(nbins);
  const double i = floor((theta - a) / (b - a) * N);

  return histogram.at(static_cast<unsigned>(i));
}

unsigned long long Histogram::total() const {
  // return std::accumulate(histogram.begin(), histogram.end(),
  //                        0, std::plus<unsigned long>());
  unsigned long long count = 0;

  for (unsigned k = 0; k < histogram.size(); k++)
    count += histogram[k];

  return count;
}

double Histogram::error() const {
  // A composition of std::transform and std::max_element might be
  // more elegant here.
  const double tot = static_cast<double>(total());
  double err = -INFINITY;               // Error in infty-norm.
  for (unsigned k = 0; k < nbins; k++) {
    const double freq = static_cast<double>(histogram[k]) / tot;
    err = std::max(err, fabs(freq - reference[k]));
  }
  return err;

#if 0
  double err = 0.0;
  for (unsigned k = 0; k < nbins; k++)
    err += fabs(histogram[k] / tot - reference[k]);
  return err / static_cast<double>(nbins);
#endif
}

std::ostream& operator<<(std::ostream& stream, const Histogram& h) {
  const auto tot = h.total();

  stream << std::setprecision(16) << tot << " ";
  for (unsigned k = 0; k < h.nbins; k++)
    stream << h.histogram[k] << " ";

  return stream << std::endl;
}

static inline double angle(unsigned k, unsigned nbins) {
  return 2.0 * M_PI / (double) nbins * (double) k - M_PI;
}

void Histogram::plot(Plotter& plotter) {
  const double tot = static_cast<double>(total());
  const double order_error = log10(error());

  std::ostringstream cmd;
  cmd << "set title 'O(Error) = " << order_error << "'\n"
      << "set xrange [-pi:pi]\n"
      << "set yrange [" << 0.0 << ":" << (maximum + minimum) << "]\n"
      << "plot '-' with linespoints pointtype 5"
      << " linecolor rgbcolor 'black' title 'Reference',"
      << " '-' with linespoints pointtype 7"
      << " linecolor rgbcolor 'blue' title 'Simulation'\n";
  for (unsigned k = 0; k < nbins; k++)
    cmd << angle(k, nbins) << " " << reference[k] << "\n";
  cmd << "\ne\n";
  for (unsigned k = 0; k < nbins; k++)
    cmd << angle(k, nbins) << " "
        << static_cast<double>(histogram[k]) / tot << "\n";
  cmd << "\ne\n";
  plotter.send(cmd.str());
}
