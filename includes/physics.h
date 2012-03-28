#ifndef PHYSICS_H
#define PHYSICS_H 1

double
magnetization(lattice_site * lattice, settings conf, gsl_vector * magnet);

double
local_energy (lattice_site * lattice, settings conf, int loc_id);

double
new_local_energy (lattice_site * lattice, settings conf, int loc_id, gsl_vector * new);

double
total_energy(lattice_site * lattice, settings conf);

#endif /* !PHYSICS_H */
