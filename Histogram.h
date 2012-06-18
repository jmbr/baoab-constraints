#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <gsl/gsl_integration.h>

#include <vector>

#include "Plotter.h"

class Histogram {
  gsl_integration_workspace* w;
  size_t nbins;
  std::vector<unsigned long> histogram;
  std::vector<double> reference;
  double temperature, K;
  simulator::Plotter plotter;
  double minimum, maximum;
  
 public:
  Histogram(size_t nbins_, double temperature_, double K_)
      : nbins(nbins_), histogram(nbins), reference(nbins),
        temperature(temperature_), K(K_) {
    w = gsl_integration_workspace_alloc(10000);
    compute();
  }

  ~Histogram() {
    gsl_integration_workspace_free(w);
  }

  unsigned long& operator[](double theta);

  double error() const;

  unsigned long long total() const;

  void plot() const;

  friend std::ostream& operator<<(std::ostream& stream, const Histogram& h);

 private:
  void compute();
};

#endif
