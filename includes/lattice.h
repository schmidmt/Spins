#ifndef LATTICE_H
#define LATTICE_H 1

#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>

gsl_vector ** allocate_lattice (int side_length, int spacedims, int spindims);

int free_lattice ( gsl_vector ** lattice , int sidelength, int spacedims );

void
print_lattice (gsl_vector ** lattice, int sidelength, int spacedims , int spindims);

void
set_homogenious_spins(gsl_vector ** lattice, int sidelength, int spacedims, int spindims);

void
set_checkerboard_spins(gsl_vector ** lattice, int sidelength, int spacedims, int spindims);

void
randomize_spins(gsl_vector ** lattice, int sidelength, int spacedims, int spindims, gsl_rng * r);

int
neighbor(int sidelength, int spacedims, int * loc, int * neigh , int num);

int *
num_to_location(int sidelength, int spacedims, int num, int * location);

int location_to_num (int sidelength, int spacedims, int * location);


#endif /* !LATTICE_H */
