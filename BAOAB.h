#ifndef BAOAB_H
#define BAOAB_H

#include <iostream>
#include <exception>

#include <gsl/gsl_rng.h>

#include <armadillo>

using namespace arma;

const double tol = 1e-14;

const unsigned nparticles = 3;                // # of particles
const unsigned nconstraints = nparticles - 1; // # of constraints

class Plotter;

class BAOAB_did_not_converge : public std::exception {
 public:
  const char* what() const throw() {
    return "Maximum number of iterations exceeded.";
  }
};

class BAOAB {
 public:
  vec::fixed<2 * nparticles> q, p;
  double friction;
  double temperature;
  double dt;

 private:
  double c1, c3;
  vec::fixed<2 * nparticles> f;

  gsl_rng* rng;

 public:
  BAOAB() {}

  BAOAB(double friction,
        double temperature,
        double dt,
        unsigned long seed);

  ~BAOAB();

  BAOAB& operator=(const BAOAB& other);

  void advance();

  double angle() const;

  double end_to_end_distance() const;

  void plot(Plotter& plotter);

 private:
  static vec::fixed<nconstraints> g(const vec& r);
  static mat::fixed<nconstraints, 2 * nparticles> G(const vec& r);

  void rattle(double h, unsigned max_iters = 1e7);
  void A();
  void compute_force();
  void B();
  void O();
};

#endif
