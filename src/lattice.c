#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <lattice.h>
#include <mikemath.h>
#include <stdlib.h>

gsl_vector **
allocate_lattice (int side_length, int spacedims, int spindims)
{
  int i;
  gsl_vector ** lattice = (gsl_vector **) malloc(intpow(side_length,spacedims)*sizeof(gsl_vector *));
  for ( i = 0 ; i < intpow(side_length,spacedims) ; i++ )
  {
    lattice[i] =  (gsl_vector *) gsl_vector_alloc (spindims);
  }
  return(lattice);
}


int
free_lattice ( gsl_vector ** lattice , int sidelength, int spacedims )
{
  int i;
  for ( i = 0 ; i < intpow(sidelength,spacedims) ; i++ )
  {
    gsl_vector_free(lattice[i]);
  }
  free(lattice);
  lattice = NULL;
  return 0;
}

void
print_lattice (gsl_vector ** lattice, int sidelength, int spacedims , int spindims)
{
  int i,j;
  int * location;
  for(i = 0 ; i < intpow(sidelength,spacedims) ; i++ )
  {
    location = num_to_location(sidelength,spacedims,i);
    printf("lattice[%d](%d) :  (",i,location_to_num(sidelength,spacedims,location));
    for(j = 0 ; j < spacedims ; j++)
    {
      printf(" %d",location[j]);
    }
    printf(" ) -> ( ");
    for(j = 0 ; j < spindims ; j++)
      printf("%e ",gsl_vector_get (lattice[i], j));
    printf(")\n");
  }
}

void
randomize_spins(gsl_vector ** lattice, int sidelength, int spacedims, int spindims)
{
  int i,j;
  const gsl_rng_type * T;
  gsl_rng * r;
  double sqrsum;
  
  gsl_rng_env_setup();
     
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for(i = 0 ; i < intpow(sidelength,spacedims) ; i++)
  {
    for(j = 0 ; j < spindims ; j++)
    {
      gsl_vector_set(lattice[i],j,2*(gsl_rng_uniform (r)-0.5));
    }
    sqrsum = 0;
    for(j = 0 ; j < spindims ; j++)
    {
      sqrsum += gsl_pow_2(gsl_vector_get(lattice[i],j));
    }
    gsl_vector_scale(lattice[i],1.0/sqrt(sqrsum));
  }

  gsl_rng_free (r);
}

gsl_vector *
neighbor(gsl_vector ** lattice, int sidelength, int spacedims)
{
  gsl_vector * neighbor = lattice[1];
  return(neighbor);
}

int *
num_to_location(int sidelength, int spacedims, int num)
{
  int * location = (int *) malloc(spacedims*sizeof(int));
  div_t asdf;
  asdf = div(num,sidelength);
  location[0] = asdf.quot;
  location[1] = asdf.rem;

  return(location);
}

int location_to_num (int sidelength, int spacedims, int * location)
{
  int num = 0;
  int i;
  for(i = 0 ; i < spacedims ; i++)
  {
    num += location[i]*intpow(sidelength,spacedims-1-i);
  }

  return(num);
}

