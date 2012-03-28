#ifndef METROPOLIS_H
#define METROPOLIS_H 1

#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <common.h>

int
mupdate_step( lattice_site * lattice, settings conf, double beta);

int
mupdatebatch(lattice_site * lattice, settings conf, double beta, datapoint * data);

#endif /* !METROPOLIS_H */
