#include <cstdio>
#include <cstdlib>

#include <getopt.h>

#include <iostream>

#include "Options.h"

static void print_help_message() {
  std::cout << "Allowed options:\n"
            << " --help                Display help information\n"
            << " --temperature arg     Temperature\n"
            << " --friction arg        Friction coefficient\n"
            << " --dt arg              Time step length\n"
            << " --min-dt arg          Minimum time step length\n"
            << " --max-dt arg          Minimum time step length\n"
            << " --number arg          Number of experiments\n"
            << " --time arg            Instant when the simulation stops\n"
            << " --seed arg            Seed for the random number generator\n"
            << " --bins arg            Number of bins in histogram\n"
            << " --plot                Plot results"
            << std::endl;
}

enum options {
  opt_help = 0,
  opt_temperature,
  opt_friction,
  opt_dt,
  opt_min_dt,
  opt_max_dt,
  opt_number,
  opt_time,
  opt_seed,
  opt_bins,
  opt_plot
};

int Options::parse(int argc, char* argv[]) {
  int c;

  while (true) {
    int option_index = 0;
    static struct option long_options[] = {
      {"help", no_argument, 0, opt_help},
      {"temperature", required_argument, 0, opt_temperature},
      {"friction", required_argument, 0, opt_friction},
      {"dt", required_argument, 0, opt_dt},
      {"min-dt", required_argument, 0, opt_min_dt},
      {"max-dt", required_argument, 0, opt_max_dt},
      {"number", required_argument, 0, opt_number},
      {"time", required_argument, 0, opt_time},
      {"seed", required_argument, 0, opt_seed},
      {"bins", required_argument, 0, opt_bins},
      {"plot", no_argument, 0, opt_plot},
      {0, 0, 0, 0}
    };

    c = getopt_long_only(argc, argv, "",
                         long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
      case opt_help:
        print_help_message();
        return -1;
        break;
      case opt_temperature:
        temperature = atof(optarg);
        break;
      case opt_friction:
        friction = atof(optarg);
        break;
      case opt_dt:
        min_dt = max_dt = atof(optarg);
        num_experiments = 1;
        break;
      case opt_min_dt:
        min_dt = atof(optarg);
        break;
      case opt_max_dt:
        max_dt = atof(optarg);
        break;
      case opt_number:
        num_experiments = atoi(optarg);
        break;
      case opt_time:
        total_time = atof(optarg);
        break;
      case opt_seed:
        random_seed = strtoul(optarg, 0, 10);
        break;
      case opt_bins:
        nbins = atoi(optarg);
        break;
      case opt_plot:
        plot = true;
        break;
      default:
        print_help_message();
        return -1;
        break;
    }
  }

  if (nbins <= 0 || total_time < 0
      || num_experiments < 1 || min_dt < 0
      || friction < 0 || temperature < 0)
  {
    std::cerr << "Not enough parameters specified. "
              << "Be sure to set at least temperature, dt, and time."
              << std::endl;
    print_help_message();
    return -1;
  }

  return 0;
}
