#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <lattice.h>
#include <mikemath.h>
#include <stdlib.h>
#include <string.h>

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
  int * location = (int *) malloc(spacedims*sizeof(int));
  for(i = 0 ; i < intpow(sidelength,spacedims) ; i++ )
  {
    location = num_to_location(sidelength,spacedims,i,location);
    printf("lattice[%d]:  (",i);
    for(j = 0 ; j < spacedims ; j++)
    {
      printf(" %d",location[j]);
    }
    printf(" ) -> ( ");
    for(j = 0 ; j < spindims ; j++)
      printf("%e ",gsl_vector_get (lattice[i], j));
    printf(")\n");
  }
  free(location);
}

void
set_homogenious_spins(gsl_vector ** lattice, int sidelength, int spacedims, int spindims)
{
  int i;
  for(i = 0 ; i < intpow(sidelength,spacedims) ; i++)
  {
    gsl_vector_set_basis(lattice[i],0);
  }
}

void
randomize_spins(gsl_vector ** lattice, int sidelength, int spacedims, int spindims, gsl_rng * r)
{
  int i,j;
  double sqrsum,random_num;
  
  for(i = 0 ; i < intpow(sidelength,spacedims) ; i++)
  {
    sqrsum = 0;
    for(j = 0 ; j < spindims ; j++)
    {
      random_num = 2*(gsl_rng_uniform(r)-0.5);
      gsl_vector_set(lattice[i],j,random_num);
      sqrsum += gsl_pow_2(random_num);
    }
    gsl_vector_scale(lattice[i],1.0/sqrt(sqrsum));
  }
}

int 
neighbor(int sidelength, int spacedims, int * loc, int * neigh , int num)
{
  //int * neigh = (int *) malloc(spacedims*sizeof(int));
  int dir = 0, ind = 0;
  if((num & 1) == 0)
    dir = -1;
  else
    dir = 1;
  
  ind = div(num,2).quot;
  if(ind >= spacedims)
  {
    fprintf(stderr,"Trying to find index out of range...\n");
    exit(EXIT_FAILURE);
  }

  memcpy(neigh,loc,spacedims*sizeof(int));

  neigh[ind] += dir;
  if(neigh[ind] < 0)
    neigh[ind] += sidelength;
  if(neigh[ind] >= sidelength)
    neigh[ind] += -sidelength;

  return(0);
}

int *
num_to_location(int sidelength, int spacedims, int num, int * location)
{
  int i,power;
  //int * location = (int *) malloc(spacedims*sizeof(int));
  div_t asdf;
  for(i = 0 ; i < spacedims ; i++)
  {
    power = intpow(sidelength,spacedims-1-i);
    asdf = div(num,power);
    location[i] = asdf.quot;
    num = asdf.rem;
  }

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
