#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <lattice.h>
#include <common.h>
#include <math.h>

double
magnetization(gsl_vector ** lattice, settings conf, gsl_vector * mag_vector)
{
  int i;
  gsl_vector_set_zero(mag_vector);
  double result = 0;
  for(i = 0 ; i < conf.elements ; i++)
  {
    gsl_vector_add(mag_vector,lattice[i]);
  }
  gsl_blas_ddot(mag_vector,mag_vector,&result);
  /*
  for(i = 0 ; i < conf.elements ; i++)
  {
    result += gsl_vector_get(lattice[i],0);
  }
  */
  return(sqrt(result));
}


double
local_energy (gsl_vector ** lattice, settings conf, int * loc)
{
  double energy = 0,result;
  int i,primary,secondary;
  int * neigh    = (int *) malloc(conf.spacedims*sizeof(int));
  primary = location_to_num(conf,loc);

  for(i = 0 ; i < 2*(conf.spacedims) ; i++)
  {
    neighbor(conf,loc,neigh,i);
    secondary = location_to_num(conf,neigh);
    gsl_blas_ddot(lattice[primary],lattice[secondary],&result);
    energy -= result;
  }
  free(neigh);
  //return(energy/(intpow(sidelength,spacedims)));
  return(energy);
}


double
new_local_energy (gsl_vector ** lattice, settings conf, int * loc, gsl_vector * new)
{
  double energy = 0,result;
  int i,secondary;
  int * neigh    = (int *) malloc(conf.spacedims*sizeof(int));
  
  for(i = 0 ; i < 2*(conf.spacedims) ; i++)
  {
    neighbor(conf,loc,neigh,i);
    secondary = location_to_num(conf,neigh);
    gsl_blas_ddot(new,lattice[secondary],&result);
    energy -= result;
  }
  free(neigh);
  //return(energy/((intpow(sidelength,spacedims))));
  return(energy);
}

double
total_energy(gsl_vector ** lattice, settings conf)
{
  int i;
  double energy = 0;
  int * loc = (int *) malloc(conf.spacedims*sizeof(int));
  for(i = 0 ; i < conf.elements ; i++)
  { 
    num_to_location(conf, i, loc);
    energy += local_energy(lattice,conf,loc)/2;
  }
  free(loc);
  return(energy);
}
