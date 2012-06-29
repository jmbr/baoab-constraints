#include <cstdlib>

#include <memory>
#include <limits>
#include <iostream>

#include <gsl/gsl_rng.h>

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
  vector<Experiment*> experiments(o.num_experiments);
  vector<double> dts = linspace<double>(o.min_dt, o.max_dt, o.num_experiments);

  gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng, o.random_seed);

  for (unsigned k = 0; k < o.num_experiments; k++) {
    const double dt = dts[k];
    experiments[k] = new Experiment(o.friction,
                                    o.temperature,
                                    dt,
                                    o.total_time,
                                    gsl_rng_get(rng),
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

  gsl_rng_free(rng);
  return EXIT_SUCCESS;
}
