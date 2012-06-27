#ifndef BAOAB_H
#define BAOAB_H

#include <iostream>
#include <exception>

#include <boost/random.hpp>

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

  boost::mt19937 rng;
  boost::random::normal_distribution<double> normal;

 public:
  BAOAB() {}

  BAOAB(double friction,
        double temperature,
        double dt,
        unsigned seed);

  BAOAB& operator=(const BAOAB& other);

  double angle() const;

  void computeForce();

  void operator()();

  void plot(Plotter& plotter);

 private:
  vec::fixed<nconstraints> g(const vec& r);
  mat::fixed<nconstraints, 2 * nparticles> G(const vec& r);

  void rattle(double h, unsigned max_iters = 1e7);
  void A();
  void B();
  void O();
};

#endif
