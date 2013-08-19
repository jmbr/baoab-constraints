#include <cassert>

#include "Force.h"

const Vector& Force::operator()(const Vector& q, double* potential) {
  force.zeros();

  force(0) =  0.0;
  force(1) = -1.0;

  pot = q(1);
  if (potential)
    *potential = pot;
  
  return force;
}
