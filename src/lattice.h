#ifndef LATTICE_H
#define LATTICE_H 1

gsl_vector **                                                                                                                                                                                  
allocate_lattice (int side_length, int spacedims, int spindims);

int
free_lattice(gsl_vector ** lattice, int sidelength, int spacedims);

void
print_lattice (gsl_vector ** lattice, int sidelength, int spacedims , int spindims);

int *
num_to_location (int sidelength, int spacedims, int num);

#endif /* !LATTICE_H */
