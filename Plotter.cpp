#include "Plotter.h"

namespace simulator {

int Plotter::open(bool persist) {
  std::string command = gnuplot_exe;
  if (persist)
    command += " -persist";
  fp = popen(command.c_str(), "w");
  return (fp != nullptr) ? 0 : -1;     // XXX Replace this by an exception.
}

void Plotter::send(std::string command) {
  // std::cout << command;
  // if (fp != nullptr)
  std::fprintf(fp, "%s", command.c_str());
  std::fflush(fp);
}

void Plotter::sendline(std::string line) {
  send(line);
  send("\n");
}

void Plotter::flush() {
  // if (fp != nullptr)
  std::fflush(fp);
}

} // namespace simulator
