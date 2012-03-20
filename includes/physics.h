#ifndef PHYSICS_H
#define PHYSICS_H 1

double
magnetization(gsl_vector ** lattice, settings conf, gsl_vector * magnet);

double
magnetization2(gsl_vector ** lattice, settings conf, gsl_vector * mag_vector);

double
local_energy (gsl_vector ** lattice, settings conf, int * loc);

double
new_local_energy (gsl_vector ** lattice, settings conf, int * loc, gsl_vector * new);

double
total_energy(gsl_vector ** lattice, settings conf);

double
total_energy2(gsl_vector ** lattice, settings conf);

#endif /* !PHYSICS_H */
