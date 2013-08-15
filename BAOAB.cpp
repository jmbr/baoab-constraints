#include <cassert>

#include <iostream>
#include <sstream>

#include "BAOAB.h"
#include "Plotter.h"
#include "abort_unless.h"

static const mat::fixed<2 * nparticles, 2 * nparticles> id =
    eye<mat>(2 * nparticles, 2 * nparticles);

BAOAB::BAOAB(double friction_, double temperature_,
             double dt_, unsigned long seed)
    : friction(friction_),
      temperature(temperature_),
      dt(dt_),
      c1(exp(-friction * dt)),
      c3(sqrt(temperature * (1.0 - c1 * c1))),
      random_seed(seed),
      rng(random_seed) {
  // Set up initial conditions (must satisfy the constraints).
  q.zeros();
  p.zeros();
  for (unsigned k = 0; k < nparticles; k++)
    q(2*k) = double(k);

  abort_unless(norm(g(q), "inf") < tol
               && norm(G(q) * p, "inf") < tol);

  f = force(q, &pot);
}

vec::fixed<nconstraints> BAOAB::g(const Vector& r) {
  vec::fixed<nconstraints> v;

  for (unsigned i = 0; i < nconstraints; i++) {
    const unsigned j = i + 1;

    vec::fixed<2> u;
    u(0) = r(2*i+0) - r(2*j+0);
    u(1) = r(2*i+1) - r(2*j+1);

    v(i) = (dot(u, u) - 1.0) / 2.0;
  }

  return v;
}

mat::fixed<nconstraints, 2 * nparticles> BAOAB::G(const Vector& r) {
  mat::fixed<nconstraints, 2 * nparticles> m;
  m.zeros();

  for (unsigned i = 0; i < nconstraints; i++) {
    const unsigned j = i + 1;

    vec::fixed<2> u;
    u(0) = r(2*i+0) - r(2*j+0);
    u(1) = r(2*i+1) - r(2*j+1);

    m(i, 2*i+0) =  u(0);
    m(i, 2*i+1) =  u(1);
    m(i, 2*i+2) = -u(0);
    m(i, 2*i+3) = -u(1);
  }

  return m;
}

void BAOAB::advance() {
  B();
  A();
  O();
  A();
  f = force(q, &pot);
  B();
}

void BAOAB::project(const Vector& q, Vector& v) {
  // Project v onto the tangent space of the manifold { x | g(x) = 0 }
  // at the point q.

  const mat::fixed<nconstraints, 2 * nparticles> Gq = G(q);
  const auto GqT = trans(Gq);
  const mat::fixed<nconstraints, nconstraints> inv_GGT = inv(symmatu(Gq * GqT));

  v -= GqT * inv_GGT * Gq * v;
}

void BAOAB::B() {
  p += dt / 2.0 * f;
  project(q, p);

  abort_unless(norm(G(q) * p, "inf") < tol);
}

void BAOAB::O() {
  for (uint k = 0; k < R.size(); ++k)
    R[k] = gaussian(rng);

  project(q, R);
  p = c1 * p + c3 * R;

  abort_unless(norm(G(q) * p, "inf") < tol);
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

    const Vector r = q - trans(Gqprev) * lambda_r;
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
  project(q, p);

  abort_unless(norm(g(q), "inf") < tol
               && norm(G(q) * p, "inf") < tol);
}

void BAOAB::plot(Plotter& plotter) {
  const double I = 5.0;
  std::ostringstream cmd;
  cmd << "unset key\n"
      << "set title 'Potential = " << pot
      << ", end-to-end distance " << end_to_end_distance() << "'\n"
      << "set xrange [" << (q(2*3+0) - I) << ":" << (q(2*3+0) + I) << "]\n"
      << "set yrange [" << (q(2*3+1) - I) << ":" << (q(2*3+1) + I) << "]\n"
      << "set size square\n"
      << "plot '-' with linespoints "
      << "pointtype 7 pointsize 3 linewidth 2\n";
  for (unsigned k = 0; k < nparticles; k++)
    cmd << q(2 * k) << " " << q(2 * k + 1) << "\n";
  cmd << "e\n";

  plotter.send(cmd.str());
}

double BAOAB::end_to_end_distance() const {
  const unsigned i = 0, j = nparticles - 1;

  vec::fixed<2> v;
  v(0) = q(2*i+0) - q(2*j+0);
  v(1) = q(2*i+1) - q(2*j+1);

  return norm(v, 2);
}

double BAOAB::potential() const {
  return pot;
}

void BAOAB::center() {
  vec::fixed<2> c = zeros<vec>(2);

  unsigned k;
  const double n = double(nparticles);

  for (k = 0; k < nparticles; k++) {
    c(0) += q(2*k+0) / n;
    c(1) += q(2*k+1) / n;
  }

  for (k = 0; k < nparticles; k++) {
    q(2*k+0) -= c(0);
    q(2*k+1) -= c(1);
  }
}
