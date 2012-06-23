#include "Plotter.h"

int Plotter::open(bool persist) {
  if (fp != nullptr)
    return 0;

  std::string command = gnuplot_exe;
  if (persist)
    command += " -persist";
  fp = popen(command.c_str(), "w");
  return (fp != 0) ? 0 : -1;     // XXX Replace this by an exception.
}

void Plotter::close() {
  if (fp != nullptr) {
    std::fclose(fp);
    fp = nullptr;
  }
}

void Plotter::send(std::string command) {
  std::fprintf(fp, "%s", command.c_str());
  std::fflush(fp);
}

void Plotter::sendline(std::string line) {
  send(line);
  send("\n");
}

void Plotter::flush() {
  // if (fp != 0)
  std::fflush(fp);
}
