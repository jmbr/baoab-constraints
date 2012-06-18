#include "Plotter.h"

namespace simulator {

int Plotter::open(bool persist) {
  std::string command = gnuplot_exe;
  if (persist)
    command += " -persist";
  fp = popen(command.c_str(), "w");
  return (fp != 0) ? 0 : -1;     // XXX Replace this by an exception.
}

void Plotter::send(std::string command) const {
  // std::cout << command;
  // if (fp != 0)
  std::fprintf(fp, "%s", command.c_str());
  std::fflush(fp);
}

void Plotter::sendline(std::string line) const {
  send(line);
  send("\n");
}

void Plotter::flush() const {
  // if (fp != 0)
  std::fflush(fp);
}

} // namespace simulator
