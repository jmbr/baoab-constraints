#ifndef BAOAB_H
#define BAOAB_H

#include <iostream>
#include <exception>
#include <random>

#include <armadillo>

#include "Misc.h"
#include "Force.h"

using namespace arma;

const double tol = 1e-12;

class Plotter;

class BAOAB_did_not_converge : public std::exception {
 public:
  const char* what() const throw() {
    return "Maximum number of iterations exceeded.";
  }
};

class BAOAB {
 public:
  Vector q, p;
  double friction;
  double temperature;
  double dt;

 private:
  double c1, c3;

  Force force;
  Vector f;
  double pot;

  unsigned long random_seed;
  std::mt19937 rng;
  std::normal_distribution<double> gaussian;

 public:
  BAOAB() {}

  BAOAB(double friction,
        double temperature,
        double dt,
        unsigned long seed);

  ~BAOAB();

  BAOAB& operator=(const BAOAB& other);

  void advance();

  void center();

  double end_to_end_distance() const;

  double potential() const;

  void plot(Plotter& plotter);

 private:
  static vec::fixed<nconstraints> g(const Vector& r);
  static mat::fixed<nconstraints, 2 * nparticles> G(const Vector& r);
  static void project(const Vector& q, Vector& v);

  void rattle(double h, unsigned max_iters = 5e7);
  void A();
  void B();
  void O();
};

#endif
