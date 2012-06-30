#ifndef BAOAB_H
#define BAOAB_H

#include <iostream>
#include <exception>

#include <gsl/gsl_rng.h>

#include <armadillo>

using namespace arma;

const double tol = 1e-12;

const unsigned nparticles = 3;                // # of particles
const unsigned nconstraints = nparticles - 1; // # of constraints

class Plotter;

typedef vec::fixed<2 * nparticles> Vector;

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
  Vector f;

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
  static void project(const Vector& q, Vector& v);

  void rattle(double h, unsigned max_iters = 5e7);
  void A();
  void compute_force();
  void B();
  void O();
};

#endif
