#ifndef METROPOLIS_H
#define METROPOLIS_H 1

#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <common.h>

int
mupdate_step( gsl_vector ** lattice, settings conf, double beta, gsl_vector * magnet, double * energy );

int
mupdate(gsl_vector ** lattice, settings conf, double beta, gsl_vector * mag_vector, double * mag, double * mag_error, double * mag2, double * mag2_err, double * energy , double * energy_error, double * energy2, double* energy2_error);

#endif /* !METROPOLIS_H */
