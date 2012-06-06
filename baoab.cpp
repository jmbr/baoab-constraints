#include <cmath>
#include <cstdlib>

#include <sys/time.h>

#include <vector>
#include <iterator>
#include <iostream>

#include <boost/random.hpp>

#include "double2.h"
#include "reference_histogram.h"
#include "Plotter.h"
#include "baoab.h"


using namespace std;

simulator::Plotter plt1, plt2;


static double D(double theta);

static long double compute_error(const vector<long double>& histogram, size_t step);

static void plot_histogram(const vector<long double>& histogram, size_t step);
static void print_histogram(const vector<long double>& histogram);

static void plot_state(const double2& q, const double2& p);

static inline unsigned get_idx(double theta) {
  const double a = -M_PI, b = M_PI;
  const double N = static_cast<double>(nbins);
  const double i = floor((theta - a) / (b - a) * N);
  return static_cast<unsigned>(i);
}


int main(int argc, char* argv[]) {
  if (argc < 6) {
    cerr << "Usage " << argv[0] << " temperature friction dt steps seed\n";
    return EXIT_FAILURE;
  }

  // Parse command line arguments.
  const double temperature = strtod(argv[1], 0);
  double friction = strtod(argv[2], 0);
  double dt = strtod(argv[3], 0);
  double steps = strtod(argv[4], 0);
  const size_t total_steps = static_cast<size_t>(steps);
  const unsigned random_seed = atoi(argv[5]);
  
  // Set up simulation.
  BAOAB baoab(friction, temperature, dt, random_seed);

  vector<long double> histogram(nbins);
  // vector<long double> bin_avg(nbins);
  // long double sqrtD_avg = 0.0;
  
  // Run simulation.
  for (size_t step = 1; step <= total_steps; step++) {
    baoab();

    const double theta = baoab.angle();
    const unsigned idx = get_idx(theta);
    ++histogram[idx];
    bin_avg[idx] += D(theta);
    sqrtD_avg += D(theta);

    if (step % static_cast<size_t>(1e5) == 0) {
      plot_histogram(histogram, step);

      print_histogram(histogram);

      // cout << (static_cast<double>(step) * dt) << " "
      //      << compute_error(histogram, step) << endl;
    }
  }

  cout << (static_cast<double>(total_steps) * dt) << " "
       << compute_error(histogram, total_steps) << endl;

  return EXIT_SUCCESS;
}


long double compute_error(const vector<long double>& histogram, size_t step) {
  // XXX Replace these lines by a call to std::accumulate()
  long double error = 0.0;

  for (size_t k = 0; k < nbins; k++)
    error += fabsl(histogram[k] / ((double) step) - reference_histogram[k]);

  return error / static_cast<long double>(nbins);
}


double D(double theta) {
  double s, c;
  sincos(theta, &s, &c);
  return sqrt(K) * c * c + s * s / sqrt(K);
}


void plot_state(const double2& q, const double2& p) {
  plt1.sendline("unset key");
  plt1.sendline("set samples 1000");
  plt1.sendline("set xrange [-2:2]");
  plt1.sendline("set yrange [-2:2]");
  plt1.sendline("set size square");
  plt1.sendline("f(x) = sqrt(1.0 - x**2)");
  plt1.send("\nplot f(x) lt -1, -f(x) lt -1, '-' with vectors filled head,");
  plt1.send(" '-' with points pointtype 7 pointsize 3,");
  plt1.sendline(" '-' with vectors filled head linetype 0");
  ostringstream os;
  os << q << " " << p << "\ne\n";
  os << q << "\ne\n";
  os << q << " " << (p / p.norm()) << "\ne\n";
  plt1.send(os.str());
}


void plot_histogram(const vector<long double>& histogram, size_t step) {
  char title[4096];
  sprintf(title, "set title 'O(Error) %Lg, steps 10^%g'\n",
          log10l(compute_error(histogram, step)), log10(step));
  plt2.send(title);
  plt2.send("plot '-' with linespoints pointtype 5 "
            "linecolor rgbcolor 'black' title 'Reference', ");
  plt2.send(" '-' with linespoints pointtype 7 "
            "linecolor rgbcolor 'blue' title 'Simulation'\n");
  for (unsigned k = 0; k < nbins; k++) {
    char buf[4096];
    sprintf(buf, "%u %Lg\n", k, reference_histogram[k]);
    plt2.send(buf);
  }
  plt2.sendline("e\n");
  for (unsigned k = 0; k < nbins; k++) {
    char buf[4096];
    sprintf(buf, "%u %Lg\n", k, histogram[k]/((double) step)) ;
    plt2.send(buf);
  }
  plt2.sendline("e\n");
}

void print_histogram(const vector<long double>& histogram) {
  ostream_iterator<double> o_i(cout, " ");
  copy(histogram.begin(), histogram.end(), o_i);
  cout << endl;
}
