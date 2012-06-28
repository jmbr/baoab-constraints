#include <iostream>
#include <sstream>

#include <gsl/gsl_randist.h>

#include "BAOAB.h"
#include "Plotter.h"
#include "abort_unless.h"

BAOAB::BAOAB(double friction_, double temperature_,
             double dt_, unsigned long seed)
    : friction(friction_),
      temperature(temperature_),
      dt(dt_),
      c1(exp(-friction * dt)),
      c3(sqrt(temperature * (1.0 - c1 * c1))) {
  rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rng, seed);

  // Set up initial conditions (must satisfy the constraints).
  q(0) = 1.0; q(1) = 0.0;
  q(2) = 0.0; q(3) = 0.0;
  q(4) =-1.0; q(5) = 0.0;

  abort_unless(norm(g(q), "inf") < tol
               && norm(G(q) * p, "inf") < tol);

  compute_force();
}

BAOAB::~BAOAB() {
  gsl_rng_free(rng);
}

BAOAB& BAOAB::operator=(const BAOAB& other) {
  if (this != &other) {
    q = other.q;
    p = other.p;
    friction = other.friction;
    temperature = other.temperature;
    dt = other.dt;
    c1 = exp(-friction * dt);
    c3 = sqrt(temperature * (1.0 - c1 * c1));
    f = other.f;
    gsl_rng_memcpy(rng, other.rng);
  }

  return *this;
}

void BAOAB::compute_force() {
  // If you change sth here, remember to change the histogram routine.
  const double d = 1.0, a = 2.0, r0 = sqrt(2.0);

  vec::fixed<2> v;
  v(0) = q(0) - q(4);
  v(1) = q(1) - q(5);

  const double r = norm(v, 2);
  const double w = exp(-a * (r - r0));
  const double morse_factor = 2.0 * a * d * (w - w * w) / r;

  f(0) = -morse_factor * v(0);
  f(1) = -morse_factor * v(1);
  f(2) =  0.0;
  f(3) =  0.0;
  f(4) =  morse_factor * v(0);
  f(5) =  morse_factor * v(1);

  // std::cout << __func__ << ": " << trans(f) << std::endl;
}

vec::fixed<nconstraints> BAOAB::g(const vec& r) {
  vec::fixed<nconstraints> v;

  for (unsigned i = 0; i < nconstraints; i++) {
    vec::fixed<2> u;
    u(0) = r(2 * i) - r(2 * i + 2);
    u(1) = r(2 * i + 1) - r(2 * i + 3);
    v(i) = (dot(u, u) - 1.0) / 2.0;
  }

  return v;
}

mat::fixed<nconstraints, 2 * nparticles> BAOAB::G(const vec& r) {
  mat::fixed<nconstraints, 2 * nparticles> m;
  m.zeros();

  for (unsigned i = 0; i < nconstraints; i++) {
    vec::fixed<2> u;

    u(0) = r(2 * i) - r(2 * i + 2);
    u(1) = r(2 * i + 1) - r(2 * i + 3);

    m(i, 2 * i)     =  u(0);
    m(i, 2 * i + 1) =  u(1);
    m(i, 2 * i + 2) = -u(0);
    m(i, 2 * i + 3) = -u(1);
  }

  return m;
}

void BAOAB::advance() {
  B();
  A();
  O();
  A();
  compute_force();
  B();
}

void BAOAB::B() {
  const mat::fixed<nconstraints, 2 * nparticles> Gq = G(q);
  const vec::fixed<nconstraints> b = Gq * (2.0 / dt * p + f);
  const vec::fixed<nconstraints> mu = solve(Gq * trans(Gq), b);

  p += dt / 2.0 * (f - trans(Gq) * mu);

  // std::cout << __func__ << ": p = " << trans(p);

  abort_unless(norm(g(q), "inf") < tol
               && norm(G(q) * p, "inf") < tol);
}

void BAOAB::O() {
  vec::fixed<2 * nparticles> R;
  for (vec::iterator r = R.begin(); r != R.end(); ++r)
    *r = gsl_ran_gaussian(rng, 1.0);

  const mat::fixed<nconstraints, 2 * nparticles> Gq = G(q);
  const vec::fixed<nconstraints> b = friction * c3 / (1.0 - c1) * G(q) * R;
  const vec::fixed<nconstraints> mu = solve(Gq * trans(Gq), b);

  p = -(1.0 - c1) / friction * trans(Gq) * mu + c1 * p + c3 * R;

  // std::cout << __func__ << ": p = " << trans(p);

  abort_unless(norm(g(q), "inf") < tol
               && norm(G(q) * p, "inf") < tol);
}

void BAOAB::A() {
  const double alpha = 1.0 / (2.0 - cbrt(2.0));
  const double beta = 1.0 - 2.0 * alpha;
  const double h = dt / 2.0;
  rattle(alpha * h);
  rattle(beta * h);
  rattle(alpha * h);
}

void BAOAB::rattle(double h, unsigned max_iters) {
  // Declare auxiliary constants.
  const mat::fixed<nconstraints, 2 * nparticles> Gqprev = G(q);

  // Deal with constraint on the configuration manifold.
  q += h * p;

  vec::fixed<nconstraints> lambda_r;
  lambda_r.zeros();

  // Solve using Newton's method.
  for (size_t k = 1; k <= max_iters; k++) {
    if (k == max_iters)
      throw BAOAB_did_not_converge();

    const vec::fixed<2 * nparticles> r = q - trans(Gqprev) * lambda_r;
    const vec::fixed<nconstraints> phi = g(r);
    const mat::fixed<nconstraints, nconstraints> dphi_dl = -G(r) * trans(Gqprev);
    const vec::fixed<nconstraints> update = solve(dphi_dl, phi);

    if (norm(phi, "inf") < tol && norm(update, "inf") < tol)
      break;

    lambda_r -= update;
  }

  q -= trans(Gqprev) * lambda_r;
  p -= trans(Gqprev) * lambda_r / h;

  // Deal with the constraint on the tangent space.
  const mat::fixed<nconstraints, 2 * nparticles> Gq = G(q);
  const vec::fixed<nconstraints> lambda_v = solve(Gq * trans(Gq), Gq * p);
  p -= trans(Gq) * lambda_v;

  abort_unless(norm(g(q), "inf") < tol
               && norm(G(q) * p, "inf") < tol);
}

void BAOAB::plot(Plotter& plotter) {
  std::ostringstream cmd;
  cmd << "unset key\n"
      << "set samples 1000\n"
      << "set xrange [-20:20]\n"
      << "set yrange [-20:20]\n"
      // << "set xrange [" << (q(2) - 2.0) << ":" << (q(2) + 2.0) << "]\n"
      // << "set yrange [" << (q(3) - 2.0) << ":" << (q(3) + 2.0) << "]\n"
      << "set size square\n"
      << "plot '-' with linespoints "
      << "pointtype 7 pointsize 1 linewidth 2\n";
  for (unsigned k = 0; k < nparticles; k++)
    cmd << q(2 * k) << " " << q(2 * k + 1) << "\n";
  cmd << "e\n";

  plotter.send(cmd.str());
}

double BAOAB::angle() const {
  // Obtain the angle between the two atoms located at the extremes of
  // the chain.
  vec::fixed<2> v10;
  v10(0) = q(0) - q(2);
  v10(1) = q(1) - q(3);

  vec::fixed<2> v12;
  v12(0) = q(4) - q(2);
  v12(1) = q(5) - q(3);

  abort_unless(fabs(norm(v10, 2) - 1.0) < tol
               && fabs(norm(v12, 2) - 1.0) < tol);

  mat::fixed<2, 2> M;
  M(0, 0) =  v10(0);
  M(0, 1) = -v10(1);
  M(1, 0) =  v10(1);
  M(1, 1) =  v10(0);

  vec::fixed<2> u = trans(M) * v12;

  return atan2(u(1), u(0));
}
