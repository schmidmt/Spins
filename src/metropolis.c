#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <lattice.h>
#include <common.h>
#include <physics.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>

int
metropolis_update(gsl_vector ** lattice, int sidelength, int spacedims, int spindims, double beta, gsl_rng * rng, gsl_vector * magnet, double * energy )
{
  int i;
  int * loc = (int *) malloc(spacedims*sizeof(int));
  double pacc,de,dep,sqrsum,random_num,exp_factor;
  gsl_vector * newspin = gsl_vector_alloc(spindims);

  /* Create a substitute random spin direction
   * which we will compare to our present spin
   */
  sqrsum = 0;
  for(i = 0 ; i < spindims ; i++)
  {
    random_num = 2*(gsl_rng_uniform(rng)-0.5);
    gsl_vector_set(newspin,i,random_num);
    sqrsum += gsl_pow_2(random_num);
  }
  gsl_vector_scale(newspin,1.0/sqrt(sqrsum));
  
  /* Choose a random location on the lattice */
  for(i = 0 ; i < spacedims ; i++)
  {
    loc[i] = gsl_rng_uniform_int(rng,sidelength);
  }
  
  de  = local_energy(lattice,sidelength,spacedims,spindims,loc);
  dep = new_local_energy(lattice,sidelength,spacedims,spindims,loc,newspin);

  /* Calculate the probibility of acceptance of the new vector, pacc.
   * If the expontential facor is less than -10, it's basically zero,
   * this code prevents it from throwing an underflow exception.
   */
  exp_factor = -(1.0/beta)*(dep-de);
  if(exp_factor > -10 && exp_factor < 1)
    pacc = gsl_sf_exp(exp_factor);
  else if(exp_factor > 1)
    pacc = 1;
  else
    pacc = 0;

  /* Generate a uniform random number between (0,1],
   * if it's < pacc, switch the vector.
   */
  random_num = gsl_rng_uniform(rng);
  if(pacc > 1 || random_num < pacc)
  {
    //gsl_vector_sub(magnet,lattice[location_to_num(sidelength,spacedims,loc)]);
    //gsl_vector_add(magnet,newspin);
    //*energy += (dep-de);
    gsl_vector_memcpy(lattice[location_to_num(sidelength,spacedims,loc)],newspin);
  }
  free(loc);
  gsl_vector_free(newspin);
  newspin = NULL;
  return(0);
}
