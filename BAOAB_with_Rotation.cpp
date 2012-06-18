#include "BAOAB.h"

void BAOAB_with_Rotation::A() {
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
