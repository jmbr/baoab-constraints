#ifndef PLOTTER_H
#define PLOTTER_H

#include <cstdio>

#include <string>

namespace simulator {

class Plotter {
  std::string gnuplot_exe;
  std::FILE* fp;

 public:
  Plotter(std::string gnuplot_exe_ = GNUPLOT_EXECUTABLE)
      : gnuplot_exe(gnuplot_exe_), fp(nullptr) {
    open();
  }
  
  ~Plotter() {
    if (fp != nullptr)
      std::fclose(fp);
  }
  
  int open(bool persist = false);

  void send(std::string command);

  void sendline(std::string line = "");

  void flush();
};

} // namespace simulator

#endif
