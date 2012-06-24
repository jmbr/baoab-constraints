#include <cstdlib>

#include <memory>
#include <limits>
#include <iostream>

#include <boost/random.hpp>

#include "Options.h"
#include "Experiment.h"
#include "linspace.h"

using namespace std;

int main(int argc, char* argv[]) {
  Options o;
  if (o.parse(argc, argv) == -1)
    return EXIT_FAILURE;

  ///////////////////////////////////////////////////////////////////////////
  // Set up experiments.
  ///////////////////////////////////////////////////////////////////////////
  boost::mt19937 rng(o.random_seed);
  const unsigned max_seed = numeric_limits<unsigned>::max();
  boost::random::uniform_int_distribution<unsigned> seeds(0, max_seed);

  vector<Experiment*> experiments(o.num_experiments);
  vector<double> dts = linspace<double>(o.min_dt, o.max_dt, o.num_experiments);

  for (unsigned k = 0; k < o.num_experiments; k++) {
    const double dt = dts[k];
    experiments[k] = new Experiment(o.K, o.friction, o.temperature, dt,
                                    o.total_time, seeds(rng), o.nbins, o.plot);
    experiments[k]->openFiles();
  }

  ///////////////////////////////////////////////////////////////////////////
  // Run simulation.
  ///////////////////////////////////////////////////////////////////////////
  unsigned r;
#pragma omp parallel for private(r)
  for (r = 0; r < o.num_experiments; r++) {
    Experiment& experiment = *experiments[r];

    try {
      experiment.simulate();
    } catch (BAOAB_did_not_converge& e) {
      cerr << "Experiment number " << r << " was not successfuly completed: "
           << e.what() << endl;
    }
  }

#pragma omp barrier

  for (unsigned k = 0; k < experiments.size(); k++) {
    experiments[k]->closeFiles();
    delete experiments[k];
  }

  return EXIT_SUCCESS;
}
