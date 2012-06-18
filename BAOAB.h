#ifndef BAOAB_H
#define BAOAB_H

#include <iostream>
#include <exception>

#include <boost/random.hpp>

#include "double2.h"
#include "Plotter.h"

using namespace std;

const double tol = 1e-14;
// const double tol = 1e-12;

class BAOAB_did_not_converge : public std::exception {
 public:
  const char* what() const throw() {
    return "Maximum number of iterations exceeded.";
  }
};

class BAOAB {
 protected:
  const double K;
  const double friction;
  const double temperature;
  const double dt;
  const double c1, c3;
  double2 f;

  // boost::mt11213b rng;
  boost::mt19937 rng;
  boost::random::normal_distribution<double> normal;

  simulator::Plotter plotter;
  
 public:
  double2 q;
  double2 p;

  BAOAB(double K_, double friction_,
        double temperature_, double dt_, unsigned seed);

  double angle() const;

  void computeForce();
  
  void operator()();

  void plot() const;
  
 protected:
  double g(const double2& r);

  double2 G(const double2& r);
  
  virtual void A() = 0;
  void B();
  void O();
};

class BAOAB_with_RATTLE : public virtual BAOAB {
 public:
  BAOAB_with_RATTLE(double K_, double friction_, double temperature_,
                    double dt_, unsigned seed)
      : BAOAB(K_, friction_, temperature_, dt_, seed) {}
  
 protected:
  void A();

 private:
  void rattle(double h, const size_t max_iters = 1e8);
};

class BAOAB_with_DoPri : public virtual BAOAB {
 public:
  BAOAB_with_DoPri(double K_, double friction_, double temperature_,
                   double dt_, unsigned seed)
      : BAOAB(K_, friction_, temperature_, dt_, seed) {}

 protected:
  void A();
};

class BAOAB_with_Rotation : public virtual BAOAB {
 public:
  BAOAB_with_Rotation(double K_, double friction_, double temperature_,
                      double dt_, unsigned seed)
      : BAOAB(1.0, friction_, temperature_, dt_, seed) {
    assert(fabs(K_ - 1.0) < 1e-32);
  }

 protected:
  void A();
};

#endif
