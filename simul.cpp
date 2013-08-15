#include <cstdlib>

#include <memory>
#include <limits>
#include <iostream>

#include "Options.h"
#include "BAOAB.h"
#include "Average.h"
#include "Plotter.h"
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
  vector<Experiment*> experiments(o.num_experiments);
  vector<double> dts = linspace<double>(o.min_dt, o.max_dt, o.num_experiments);

  std::mt19937 rng;
  rng.seed(o.random_seed);
  std::uniform_int_distribution<unsigned long> uniform;

  for (unsigned k = 0; k < o.num_experiments; k++) {
    const double dt = dts[k];
    const unsigned long seed = uniform(rng);
    clog << "Initializing random number generator for "
         << "simulation using dt = " << dt
         << " with seed " << seed << endl;
    experiments[k] = new Experiment(o.friction,
                                    o.temperature,
                                    dt,
                                    o.equilibration_time,
                                    o.production_time,
                                    seed,
                                    o.plot);
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
