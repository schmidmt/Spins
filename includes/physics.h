#ifndef PHYSICS_H
#define PHYSICS_H 1

double
magnetization(gsl_vector ** lattice, settings conf, gsl_vector * magnet);

double
local_energy (gsl_vector ** lattice, settings conf, int * loc);

double
new_local_energy (gsl_vector ** lattice, settings conf, int * loc, gsl_vector * new);

double
total_energy(gsl_vector ** lattice, settings conf);

#endif /* !PHYSICS_H */
