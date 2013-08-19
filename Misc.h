#ifndef MISC_H
#define MISC_H

#include <armadillo>

using namespace arma;

const double K = 0.05;

const unsigned nparticles = 1;                // # of particles
const unsigned nconstraints = 1;              // # of constraints

typedef vec::fixed<2 * nparticles> Vector;

#endif
