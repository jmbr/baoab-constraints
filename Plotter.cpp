#include "Plotter.h"

int Plotter::open(bool persist) {
  if (fp != 0)
    return 0;

  std::string command = gnuplot_exe;
  if (persist)
    command += " -persist";
  fp = popen(command.c_str(), "w");
  return (fp != 0) ? 0 : -1;     // XXX Replace this by an exception.
}

void Plotter::close() {
  if (fp != 0) {
    std::fclose(fp);
    fp = 0;
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
