#ifndef LATTICE_H
#define LATTICE_H 1

#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <common.h>

gsl_vector **
allocate_lattice (settings conf);

int
free_lattice ( gsl_vector ** lattice , settings conf );

void
print_lattice (gsl_vector ** lattice, settings conf );

void
set_homogenious_spins(gsl_vector ** lattice, settings conf );

void
set_checkerboard_spins(gsl_vector ** lattice, settings conf);

void
randomize_spins(gsl_vector ** lattice, settings conf );

int 
neighbor(settings conf, int * loc, int * neigh , int num);

int *
num_to_location(settings conf, int num, int * location);

int
location_to_num(settings conf, int * location);


#endif /* !LATTICE_H */
