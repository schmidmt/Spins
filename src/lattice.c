#include <lattice.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <math.h>

gsl_vector **
allocate_lattice (int side_length, int spacedims, int spindims)
{
  int i;
  gsl_vector ** lattice = (gsl_vector **) malloc((side_length**spacedims)*sizeof(gsl_vector));
  for ( i = 0 ; i < side_length**spacedims ; i++ )
  {
    lattice[i] =  (gsl_vector *) gsl_vector_alloc (spindims);
  }
  return();
}


