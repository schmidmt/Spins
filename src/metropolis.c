#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <lattice.h>
#include <mikemath.h>
#include <physics.h>
#include <gsl/gsl_sf_exp.h>


int
metropolis_update(gsl_vector ** lattice, int sidelength, int spacedims, int spindims, double beta, gsl_rng * rng)
{
  int i;
  int * loc = (int *) malloc(spacedims*sizeof(int));
  double pacc,de,dep,sqrsum,random_num;
  gsl_vector * newspin = gsl_vector_alloc(spindims);
  sqrsum = 0;
  for(i = 0 ; i < spindims ; i++)
  {
    random_num = 2*(gsl_rng_uniform(rng)-0.5);
    gsl_vector_set(newspin,i,random_num);
    sqrsum += gsl_pow_2(random_num);
  }
  gsl_vector_scale(newspin,1.0/sqrt(sqrsum));
  
  for(i = 0 ; i < spacedims ; i++)
  {
    loc[i] = gsl_rng_uniform_int(rng,sidelength);
  }
  
  de  = local_energy(lattice,sidelength,spacedims,spindims,loc);
  dep = new_local_energy(lattice,sidelength,spacedims,spindims,loc,newspin);
  
  /*
  printf("Metropolis_Update: Lattice Point = (");
  for(i = 0 ; i < spacedims ; i++)
  {
    printf(" %d",loc[i]);
  }
  printf(" ) : de = %e  dep = %e",de,dep);
  */
  pacc = gsl_sf_exp(-beta*(dep-de));
  if(pacc > 1)
    pacc = 1;

  if(pacc == 1 || gsl_rng_uniform(rng) < pacc)
  {
    /*printf(" A\n");*/
    gsl_vector_memcpy(lattice[location_to_num(sidelength,spacedims,loc)],newspin);
  }
   
  else
  {
    /*printf("\n");*/
  }
  free(loc);
  gsl_vector_free(newspin);
  newspin = NULL;
  return(0);
}
