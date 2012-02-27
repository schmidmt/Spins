#ifndef METROPOLIS_H
#define METROPOLIS_H 1

#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>

int
metropolis_update ( gsl_vector ** lattice, int sidelength, int spacedims, int spindims, double beta, gsl_rng * rng, gsl_vector * magnet, double * energy );

#endif /* !METROPOLIS_H */
