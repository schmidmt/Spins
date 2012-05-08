#ifndef BLOCKING_H
#define BLOCKING_H 1

#include <gsl/gsl_vector.h>
#include <common.h>

int
block_4_majority(settings confin,lattice_site * lin , settings confout, \
                 lattice_site * lout);



#endif /* !BLOCKING_H */
