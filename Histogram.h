#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <vector>

class Plotter;

class Histogram {
  unsigned nbins;
  std::vector<unsigned long> histogram;
  std::vector<double> reference;
  double temperature;
  double minimum, maximum;

 public:
  Histogram() {}

  Histogram(unsigned nbins_, double temperature_);

  unsigned long& operator[](double theta);

  double error() const;

  unsigned long long total() const;

  void plot(Plotter& plotter);

  inline unsigned size() const { return nbins; }

  friend std::ostream& operator<<(std::ostream& stream, const Histogram& h);

 private:
  void compute();
};

#endif
