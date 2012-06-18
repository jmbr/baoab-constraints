#include <cstdlib>

#include <memory>
#include <limits>
#include <iostream>

#include <boost/random.hpp>

#include "Options.h"
#include "Experiment.h"
#include "linspace.h"

using namespace std;
using namespace simulator;

int main(int argc, char* argv[]) {
  Options o;
  if (o.parse(argc, argv) == -1)
    return EXIT_FAILURE;

  ///////////////////////////////////////////////////////////////////////////
  // Set up experiments.
  ///////////////////////////////////////////////////////////////////////////
  boost::mt19937 rng(o.random_seed);
  const unsigned s = numeric_limits<unsigned>::max();
  boost::random::uniform_int_distribution<unsigned> seeds(0, s);
  
  vector<unique_ptr<Experiment> > experiments(o.num_experiments);
  vector<double> dts = linspace<double>(o.min_dt, o.max_dt, o.num_experiments);
  
  for (size_t k = 0; k < experiments.size(); k++) {
    const double dt = dts[k];
    experiments[k].reset(new Experiment(o.K, o.friction, o.temperature, dt,
                                        o.total_time, seeds(rng), o.nbins,
                                        o.plot));
  }

  ///////////////////////////////////////////////////////////////////////////
  // Run simulation.
  ///////////////////////////////////////////////////////////////////////////
  size_t r;
#pragma omp parallel for private(r)
  for (r = 0; r < o.num_experiments; r++) {
    Experiment& experiment = *experiments[r].get();

    try {
      experiment.simulate();
    } catch (BAOAB_did_not_converge& e) {
      cerr << "Experiment number " << r << " was not successfuly completed: "
           << e.what() << endl;
    }
  }
#pragma omp barrier
  
  return EXIT_SUCCESS;
}
