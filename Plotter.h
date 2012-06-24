#ifndef PLOTTER_H
#define PLOTTER_H

#include <cstdio>

#include <string>

class Plotter {
  std::string gnuplot_exe;
  std::FILE* fp;

 public:
  Plotter(std::string gnuplot_exe_ = GNUPLOT_EXECUTABLE)
      : gnuplot_exe(gnuplot_exe_), fp(0) {}

  ~Plotter() {
    close();
  }

  int open(bool persist = false);

  void close();

  void send(std::string command);

  void sendline(std::string line = "");

  void flush();
};

#endif
