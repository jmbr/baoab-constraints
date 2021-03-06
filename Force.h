#ifndef FORCE_H
#define FORCE_H

#include "Misc.h"

class Force {
 public:
  Force() : pot(0.0) {}

  const Vector& operator()(const Vector& q, double* potential = 0);
  
 private:
  double pot;
  Vector force;
};

#endif
