#include <cassert>

#include "Force.h"

static const double interactions_table[8 * 7] = { // The 8 is for alignment purposes.
  0.000000, 1.000000, 1.732051, 1.000000, 2.000000, 1.732051, 1.000000, 0.000000,
  1.000000, 0.000000, 1.000000, 1.000000, 1.732051, 2.000000, 1.732051, 0.000000,
  1.732051, 1.000000, 0.000000, 1.000000, 1.000000, 1.732051, 2.000000, 0.000000,
  1.000000, 1.000000, 1.000000, 0.000000, 1.000000, 1.000000, 1.000000, 0.000000,
  2.000000, 1.732051, 1.000000, 1.000000, 0.000000, 1.000000, 1.732051, 0.000000,
  1.732051, 2.000000, 1.732051, 1.000000, 1.000000, 0.000000, 1.000000, 0.000000,
  1.000000, 1.732051, 2.000000, 1.000000, 1.732051, 1.000000, 0.000000, 0.000000,
};

static inline unsigned idx(unsigned i, unsigned j) { return 8*i+j; }

void Force::compute_force_morse(const Vector& q, unsigned i, unsigned j) {
  assert(j - i > 1);
  
  const double d = 1.0, a = 2.0;
  const double r0 = interactions_table[idx(i, j)];

  vec::fixed<2> v;
  v(0) = q(2*i+0) - q(2*j+0);
  v(1) = q(2*i+1) - q(2*j+1);

  const double r = norm(v, 2);
  const double w = exp(-a * (r - r0));
  const double morse_factor = 2.0 * a * d * (w - w * w) / r;

  pot += d * (1.0 - w) * (1.0 - w);
  force(2*i+0) -= morse_factor * v(0);
  force(2*i+1) -= morse_factor * v(1);
  force(2*j+0) += morse_factor * v(0);
  force(2*j+1) += morse_factor * v(1);
}

const Vector& Force::operator()(const Vector& q, double* potential) {
  pot = 0.0;
  force.zeros();
  
  for (unsigned i = 0; i < nparticles; i++)
    for (unsigned j = i+2; j < nparticles; j++)
      if (j - i > 1)
        compute_force_morse(q, i, j);

  if (potential)
    *potential = pot;
  
  return force;
}
