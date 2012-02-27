#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <lattice.h>
#include <common.h>
#include <math.h>

double
magnetization(gsl_vector ** lattice, int sidelength, int spacedims, int spindims, gsl_vector * magnet)
{
  int i;
  gsl_vector_set_zero(magnet);
  double result;
  for(i = 0 ; i < intpow(sidelength,spacedims) ; i++)
  {
    gsl_vector_add(magnet,lattice[i]);
  }
  gsl_blas_ddot(magnet,magnet,&result);
  return(sqrt(result));
}


double
local_energy (gsl_vector ** lattice, int sidelength, int spacedims, int spindims, int * loc)
{
  double energy = 0,result;
  int i,primary,secondary;
  int * neigh    = (int *) malloc(spacedims*sizeof(int));
  primary = location_to_num(sidelength,spacedims,loc);

  for(i = 0 ; i < 2*spacedims ; i++)
  {
    neighbor(sidelength,spacedims,loc,neigh,i);
    secondary = location_to_num(sidelength,spacedims,neigh);
    gsl_blas_ddot(lattice[primary],lattice[secondary],&result);
    energy += result;
  }
  free(neigh);
  //return(energy/(intpow(sidelength,spacedims)));
  return(energy);
}


double
new_local_energy (gsl_vector ** lattice, int sidelength, int spacedims, int spindims, int * loc, gsl_vector * new)
{
  double energy = 0,result;
  int i,secondary;
  int * neigh    = (int *) malloc(spacedims*sizeof(int));
  
  for(i = 0 ; i < 2*spacedims ; i++)
  {
    neighbor(sidelength,spacedims,loc,neigh,i);
    secondary = location_to_num(sidelength,spacedims,neigh);
    gsl_blas_ddot(new,lattice[secondary],&result);
    energy += result;
  }
  free(neigh);
  //return(energy/((intpow(sidelength,spacedims))));
  return(energy);
}

double
total_energy(gsl_vector ** lattice, int sidelength, int spacedims, int spindims )
{
  int i;
  double energy = 0;
  int * loc = (int *) malloc(spacedims*sizeof(int));
  for(i = 0 ; i < intpow(sidelength,spacedims) ; i++)
  { 
    num_to_location(sidelength, spacedims, i, loc);
    energy += local_energy(lattice,sidelength,spacedims,spindims,loc);
  }
  return(energy);
}
