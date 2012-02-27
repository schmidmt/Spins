#ifndef PHYSICS_H
#define PHYSICS_H 1

double
magnetization(gsl_vector ** lattice, int sidelength, int spacedims, int spindims, gsl_vector * magnet);

double
local_energy (gsl_vector ** lattice, int sidelength, int spacedims, int spindims, int * loc);

double
new_local_energy (gsl_vector ** lattice, int sidelength, int spacedims, int spindims, int * loc, gsl_vector * new);

double
total_energy(gsl_vector ** lattice, int sidelength, int spacedims, int spindims );

#endif /* !PHYSICS_H */
