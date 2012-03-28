#ifndef LATTICE_H
#define LATTICE_H 1

#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <common.h>

lattice_site *
allocate_lattice (settings conf);

int
free_lattice ( lattice_site * lattice , settings conf );

void
print_lattice (lattice_site * lattice, settings conf );

void
set_homogenious_spins(lattice_site * lattice, settings conf );

void
set_checkerboard_spins(lattice_site * lattice, settings conf);

int
get_neighbor_id(settings conf, int site_id, int num);

void
randomize_spins(lattice_site * lattice, settings conf );

int 
find_neighbors(settings conf, int * loc, int * neigh , int num);

inline int *
num_to_location(settings conf, int num, int * location);

inline int
location_to_num(settings conf, int * location);


#endif /* !LATTICE_H */
