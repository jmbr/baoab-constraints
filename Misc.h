#ifndef MISC_H
#define MISC_H

#include <armadillo>

using namespace arma;

const unsigned nparticles = 7;                // # of particles
const unsigned nconstraints = nparticles - 1; // # of constraints

typedef vec::fixed<2 * nparticles> Vector;

#endif
