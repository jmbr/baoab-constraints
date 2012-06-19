#include <iostream>

#include <boost/program_options.hpp>

#include "Options.h"

namespace po = boost::program_options;

int Options::parse(int argc, char* argv[]) {
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help", "Display help information")
      ("K", po::value<double>(), "Length of the ellipse's horizontal axis")
      ("temperature", po::value<double>(), "Temperature")
      ("friction", po::value<double>(), "Friction coefficient")
      ("dt", po::value<double>(), "Time step length")
      ("min-dt", po::value<double>(), "Minimum time step length")
      ("max-dt", po::value<double>(), "Minimum time step length")
      ("number", po::value<size_t>(), "Number of experiments")
      ("time", po::value<double>(), "Instant when the simulation stops")
      ("seed", po::value<unsigned>(), "Seed for the random number generator")
      ("bins", po::value<unsigned>(), "Number of bins in histogram")
      ("plot", po::value<bool>(), "Plot results")
      ;

  po::variables_map varmap;
  po::store(po::parse_command_line(argc, argv, desc), varmap);
  po::notify(varmap);

  if (varmap.count("help")) {
    std::cout << desc;
    return -1;
  }

  if (varmap.count("K"))
    K = varmap["K"].as<double>();

  if (varmap.count("temperature"))
    temperature = varmap["temperature"].as<double>();

  if (varmap.count("friction"))
    friction = varmap["friction"].as<double>();

  if (varmap.count("dt")) {
    min_dt = varmap["dt"].as<double>();
    max_dt = min_dt;
    num_experiments = 1;
  } else {
    if (varmap.count("min-dt"))
      min_dt = varmap["min-dt"].as<double>();

    if (varmap.count("max-dt"))
      max_dt = varmap["max-dt"].as<double>();

    if (varmap.count("number"))
      num_experiments = varmap["number"].as<size_t>();
  }

  if (varmap.count("time"))
    total_time = varmap["time"].as<double>();

  if (varmap.count("seed"))
    random_seed = varmap["seed"].as<unsigned>();

  if (varmap.count("bins"))
    nbins = varmap["bins"].as<unsigned>();

  if (varmap.count("plot"))
    plot = varmap["plot"].as<bool>();

  if (nbins <= 0 || total_time < 0
      || num_experiments < 1 || min_dt < 0
      || friction < 0 || temperature < 0 || K < 0)
  {
    std::cerr << desc;
    return -1;
  }

  return 0;
}
