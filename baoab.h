#ifndef BAOAB_H
#define BAOAB_H

#include <iostream>

#include "double2.h"

extern "C" {
// #include "dopri5.h"
#include "dop853.h"
}

// const double K = 3.0 / 2.0;
const double K = 5.0;
// const double K = 1.0;

const double tol = 1e-11;

using namespace std;

class BAOAB {
  const double friction;
  const double temperature;
  const double dt;
  const double c1, c3;
  double2 f;

  boost::mt11213b rng;
  // boost::mt19937 rng;
  boost::random::normal_distribution<double> normal;

 public:
  double2 q;
  double2 p;

  BAOAB(double friction_, double temperature_,
        double dt_, unsigned seed)
      : friction(friction_),
        temperature(temperature_),
        dt(dt_),
        c1(exp(-friction * dt)),
        c3(sqrt(temperature * (1.0 - c1 * c1))) {
    rng.seed(seed);

    // Set up initial conditions (must satisfy the constraints).
    // boost::random::uniform_real_distribution<double> unif(-1.0, 1.0);
    // q.x = unif(rng);
    // q.y = unif(rng);
    // q = q / q.norm();
    // p.x = -q.y;
    // p.y =  q.x;

    q.x = 0;
    q.y = 1;
    p.x = 1;
    p.y = 0;

    assert(fabs(g(q)) < tol);
    assert(fabs(dot(G(q), p)) < tol);
    
    computeForce();
  }

  inline void computeForce() {
    // f.x = -q.y * 2.0 * q.x;
    // f.y = -q.x * q.x;
    f.x =  0.0;
    f.y = -1.0;
  }
  
  void operator()() {
    B();
    A();
    O();
    A();
    computeForce();
    B();
  }

  inline double angle() const {
    return atan2(q.y, q.x);
  }
  
 private:
  static inline double g(const double2& r) {
    return K * r.x * r.x + r.y * r.y - 1.0;
  }

  static inline double2 G(const double2& r) {
    return double2(2.0 * K * r.x, 2.0 * r.y);
  }
  
  void B() {
    // std::cerr << "B" << std::endl;
  
    const double2 Gq = G(q);
    const double Gq_times_Gq_trans = dot(Gq, Gq);
    const double2 aux = 2.0 / dt * p + f;
    const double mu = dot(Gq, aux) / (Gq_times_Gq_trans);
    p += dt / 2.0 * (f - mu * Gq);

    assert(fabs(g(q)) < tol && fabs(dot(G(q), p)) < tol);
  }

  /*
  static void rattle(double2& q, double2& p, double h,
                     const size_t max_iters = 1e5) {
    // Declare auxiliary constants.
    const double2 Gqprev = G(q);
    
    // Deal with constraint on the configuration manifold.
    q += h * p;
    double lambda_r = 0.0;
    // Solve using Newton's method.
    for (size_t k = 1; k <= max_iters; k++) {
      if (k == max_iters) {
        std::cerr << "Maximum number of iterations exceeded." << std::endl;
        // std::cerr << "Orig q = " << qprev << " q = " << q << " h = " << h << std::endl;
        abort();
      }

      const double2 r = q - lambda_r * Gqprev;
      const double phi = g(r);
      const double dphi_dl = -dot(G(r), Gqprev);
      const double update = phi / dphi_dl;
      
      if (fabs(phi) < tol && fabs(update) < tol)
        break;
      
      lambda_r -= update;
    }

    q -= lambda_r * Gqprev;
    p -= lambda_r / h * Gqprev;

    // Deal with constraint on the tangent space.
    const double2 Gq = G(q);
    double lambda_v = dot(Gq, p) / dot(Gq, Gq);
    p -= lambda_v * Gq;

    assert(fabs(g(q)) < tol && fabs(dot(G(q), p)) < tol);
  }

  void A() {
    // std::cerr << "A" << endl;

    const double h = dt / 2.0;

    const double alpha = 1.0 / (2.0 - cbrt(2.0));
    const double beta = 1.0 - 2.0 * alpha;
    rattle(q, p, alpha * h);
    rattle(q, p,  beta * h);
    rattle(q, p, alpha * h);

    // std::cout << fabs(g(q)) << " " << fabs(dot(G(q), p)) << std::endl;
  }
  */

  /*
  void A() {
    // std::cerr << "A" << endl;

    const double h = dt / 2.0;

    const size_t n = 10;
    for (size_t k = 0; k < 10; k++)
      rattle(q, p, h / static_cast<double>(n));
  }
  */

  /*
  void A() {
    // std::cerr << "A" << endl;

    double theta = angle();
    const double dtheta = q.x * p.y - q.y * p.x;

    theta += dt / 2.0 * dtheta;

    double s, c;
    sincos(theta, &s, &c);

    q.x =  c;
    q.y =  s;

    p.x = -s * dtheta;
    p.y =  c * dtheta;

    assert(fabs(g(q)) < tol);
    assert(fabs(dot(G(q), p)) < tol);
  }
  */
  
  static void fellipse(unsigned n __attribute__((unused)),
                       double x __attribute__((unused)),
                       double* y,
                       double* f) {
    const double two_lambda = (K * y[2] * y[2] + y[3] * y[3]) / (K * K * y[0] * y[0] + y[1] * y[1]);
    f[0] = y[2];
    f[1] = y[3];
    f[2] = -K * y[0] * two_lambda;
    f[3] = -y[1] * two_lambda;
  }

  void A() {
    // std::cerr << "A" << std::endl;

    // if (fabs(g(q)) >= tol || fabs(dot(G(q), p)) >= tol)
    //   cerr << "g(q) = " << g(q) << " " << "G(q) p = " << dot(G(q), p) << endl;

    // assert(fabs(g(q)) < tol);
    // assert(fabs(dot(G(q), p)) < tol);

    unsigned ndgl = 4;
    unsigned licont = 2;
    double y[ndgl];
    // unsigned icont[licont], i;
    // int res;
    int iout, itoler;
    double x, xend, atoler, rtoler;

    iout = 0;                      // Don't show output for each step.
    x = 0.0;
    y[0] = q.x;
    y[1] = q.y;
    y[2] = p.x;
    y[3] = p.y;
    xend = dt / 2.0;                    // <--- Advance one half step.
    itoler = 0;
    rtoler = 1e-16;
    atoler = rtoler;
    // icont[0] = 0;
    // icont[1] = 1;

    // res =
    // dopri5 (ndgl, fellipse, x, y, xend, &rtoler, &atoler, itoler, solout, iout,
    dop853 (ndgl, fellipse, x, y, xend, &rtoler, &atoler, itoler, nullptr, iout,
            stderr, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, ndgl, nullptr, licont);

    q.x = y[0];
    q.y = y[1];
    p.x = y[2];
    p.y = y[3];

    // std::cout << fabs(g(q)) << " " << fabs(dot(G(q), p)) << std::endl;
    
    // if (fabs(g(q)) >= tol || fabs(dot(G(q), p)) >= tol)
    //   cerr << "g(q) = " << g(q) << " " << "G(q) p = " << dot(G(q), p) << endl;
    
    assert(fabs(g(q)) < tol && fabs(dot(G(q), p)) < tol);
  }

  /*
  void O() {
    // std::cerr << "O" << std::endl;
    
    const double x = q.x, y = q.y;
    const double K2 = K * K, x2 = x * x, y2 = y * y;
    const double denom = x2 * K2 + y2;
    // M = I - H, where H = G' * G / (G * G')
    const double2 Mrow1(1.0 - x2 * K2 / denom, - x * y * K / denom);
    const double2 Mrow2(-x * y * K / denom, 1.0 - y2 / denom);

    const double2 S(normal(rng), normal(rng));
    const double2 R(dot(Mrow1, S), dot(Mrow2, S));
    
    p = c1 * p + c3 * R;

    assert(fabs(g(q)) < tol && fabs(dot(G(q), p)) < tol);
  }
  */
  
  void O() {
    // std::cerr << "O" << std::endl;

    const double2 Gq = G(q);
    const double2 R = { normal(rng), normal(rng) };
    const double2 aux = c1 * p + c3 * R;
    const double mu = dot(Gq, aux) / ((1 - c1) / friction * dot(Gq, Gq));
    const double2 psi = -Gq * mu / friction;

    p = psi + c1 * (psi - p) + c3 * R;

    assert(fabs(g(q)) < tol && fabs(dot(G(q), p)) < tol);
  }
};

#endif
