#include "BAOAB.h"

void BAOAB_with_RATTLE::A() {
  const double alpha = 1.0 / (2.0 - cbrt(2.0));
  const double beta = 1.0 - 2.0 * alpha;
  const double h = dt / 2.0;
  rattle(alpha * h);
  rattle(beta * h);
  rattle(alpha * h);
}

/*
void BAOAB_with_RATTLE::A() {
  const double h = dt / 2.0;

  const size_t n = 10;
  for (size_t k = 0; k < 10; k++)
    rattle(q, p, h / static_cast<double>(n));
}
*/

void BAOAB_with_RATTLE::rattle(double h, const size_t max_iters) {
  // Declare auxiliary constants.
  const double2 Gqprev = G(q);

  // Deal with constraint on the configuration manifold.
  q += h * p;
  double lambda_r = 0.0;
  // Solve using Newton's method.
  for (size_t k = 1; k <= max_iters; k++) {
    if (k == max_iters)
      throw BAOAB_did_not_converge();

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
