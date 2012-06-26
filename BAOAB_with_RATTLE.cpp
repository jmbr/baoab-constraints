#include "BAOAB.h"

#include "abort_unless.h"

void BAOAB_with_RATTLE::A() {
  const double alpha = 1.0 / (2.0 - cbrt(2.0));
  const double beta = 1.0 - 2.0 * alpha;
  const double h = dt / 2.0;
  rattle(alpha * h);
  rattle(beta * h);
  rattle(alpha * h);
}

void BAOAB_with_RATTLE::rattle(double h, unsigned max_iters) {
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
