#include <iostream>

extern "C" {
#include "dop853.h"
// #include "dopri5.h"
}

#include "BAOAB.h"

static void fellipse(unsigned n __attribute__((unused)),
                     double x __attribute__((unused)),
                     double* y,
                     double* func,
                     double K) {
  const double K_sqr = K * K;
  const double y0_sqr = y[0] * y[0];
  const double y1_sqr = y[1] * y[1];
  const double y2_sqr = y[2] * y[2];
  const double y3_sqr = y[3] * y[3];
  const double two_lambda = (K * y2_sqr + y3_sqr) / (K_sqr * y0_sqr + y1_sqr);
  func[0] = y[2];
  func[1] = y[3];
  func[2] = -K * y[0] * two_lambda;
  func[3] = -y[1] * two_lambda;
}

void BAOAB_with_DoPri::A() {
  unsigned ndgl = 4;
  unsigned licont = 2;
  double y[ndgl];
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
  rtoler = 1e-19;
  atoler = rtoler;

  // dopri5(ndgl, fellipse, x, y, xend, &rtoler, &atoler, itoler, 0, iout,
  //        stderr, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, ndgl, 0, licont, K);
  dop853(ndgl, fellipse, x, y, xend, &rtoler, &atoler, itoler, 0, iout,
         stderr, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100000000, 0, 0, ndgl, 0, licont, K);

  q.x = y[0];
  q.y = y[1];
  p.x = y[2];
  p.y = y[3];

#if 0
  static unsigned lll = 0;
  if (lll % 1000 == 0)
    std::cerr << "log10(fabs(g(q))) = " << log10(fabs(g(q))) << std::endl;
  ++lll;
  if (fabs(g(q)) >= tol) {
    std::cerr << "g(q) = " << g(q) << std::endl;
    abort();
  }
  if (fabs(dot(G(q), p)) >= tol) {
    std::cerr << "G(q) p = " << dot(G(q), p) << std::endl;
    abort();
  }
#endif

  assert(fabs(g(q)) < tol);
  assert(fabs(dot(G(q), p)) < tol);
}
