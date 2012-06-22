#include <iostream>
#include <sstream>

#include "BAOAB.h"

BAOAB::BAOAB(double K_,
             double friction_, double temperature_,
             double dt_, unsigned seed)
    : K(K_),
      friction(friction_),
      temperature(temperature_),
      dt(dt_),
      c1(exp(-friction * dt)),
      c3(sqrt(temperature * (1.0 - c1 * c1))) {
  rng.seed(seed);

  // Set up initial conditions (must satisfy the constraints).
  q.x = 0.0;
  q.y = 1.0;
  p.x = -1.0;
  p.y = 0.0;
  // boost::random::uniform_real_distribution<double> unif0_2pi(0.0, 2.0 * M_PI);
  // const double theta = unif0_2pi(rng);
  // q.x = cos(theta) / sqrt(K);
  // q.y = sin(theta);
  // p.x = -q.y;
  // p.y = K * q.x;

  assert(fabs(g(q)) < tol);
  assert(fabs(dot(G(q), p)) < tol);

  computeForce();
}

void BAOAB::computeForce() {
  // If you change sth here, remember to change the histogram routine.

  const double d = 1.0, a = 2.0, r0 = 1.0, r = q.norm();
  const double w = exp(-a * (r - r0));
  const double morse_factor = -2.0 * a * d * (w - w * w) / r;
  f.x = morse_factor * q.x;
  f.y = morse_factor * q.y;

  // f.x = -q.y * 2.0 * q.x; f.y = -q.x * q.x;
  // f.x = -2.0 * q.x * q.y * q.y; f.y = -2.0 * q.x * q.x * q.y;
  // f.x =  0.0; f.y = -1.0;
}

void BAOAB::operator()() {
  B();
  A();
  O();
  A();
  computeForce();
  B();
}

double BAOAB::angle() const {
  return atan2(q.y, sqrt(K) * q.x);
}

double BAOAB::g(const double2& r) {
  return K * r.x * r.x + r.y * r.y - 1.0;
}

double2 BAOAB::G(const double2& r) {
  return double2(2.0 * K * r.x, 2.0 * r.y);
}

void BAOAB::B() {
  const double2 Gq = G(q);
  const double Gq_times_Gq_trans = dot(Gq, Gq);
  const double2 aux = 2.0 / dt * p + f;
  const double mu = dot(Gq, aux) / (Gq_times_Gq_trans);
  p += dt / 2.0 * (f - mu * Gq);

  assert(fabs(g(q)) < tol && fabs(dot(G(q), p)) < tol);
}

void BAOAB::O() {
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

void BAOAB::plot() const {
  ostringstream cmd;
  cmd << "unset key\n"
      << "set samples 1000\n"
      << "set xrange [-2:2]\n"
      << "set yrange [-2:2]\n"
      << "plot '-' with points pointtype 7 pointsize 3,"
      << " '-' with vectors filled head linetype 0\n"
      << q << "\ne\n"
      << q << " " << (p / p.norm()) << "\ne\n";
  plotter.send(cmd.str());
}
