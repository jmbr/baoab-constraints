#include <cassert>

#include "Average.h"

Average::Average() : count(0ULL), sum(0.0), c(0.0) {}

void Average::update(double input) {
  // Implementation of Kahan summation (see
  // http://en.wikipedia.org/wiki/Kahan_summation)
  ++count;
  long double y = input - c;
  const long double t = sum + y;
  c = (t - sum) - y;
  sum = t;
}

long double Average::value() const {
  assert(count > 0);
  return sum / static_cast<long double>(count);
}
