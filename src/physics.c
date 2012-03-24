/******************************************************************************
    Copyright 2012 Michael Schmidt (mts@colorado.edu)

    This file is part of spins.

    spins is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    spins is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with spins.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/


#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <lattice.h>
#include <common.h>
#include <math.h>

double
magnetization2(gsl_vector ** lattice, settings conf, gsl_vector * mag_vector)
{
  int i;
  gsl_vector_set_zero(mag_vector);
  double result = 0;
  double mag2 = 0;
  for(i = 0 ; i < conf.elements ; i++)
  {
    gsl_blas_ddot(lattice[i],lattice[i],&result);
    mag2 += result;
  }
  return(mag2/conf.elements);
}

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
  return(result);
  */
  /* I have no idea why this is happening! */
  return(result/intpow(conf.elements,2));
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
  return(energy/conf.elements);
}

double
total_energy2(gsl_vector ** lattice, settings conf)
{
  int i;
  double energy2 = 0;
  int * loc = (int *) malloc(conf.spacedims*sizeof(int));
  for(i = 0 ; i < conf.elements ; i++)
  { 
    num_to_location(conf, i, loc);
    energy2 += gsl_pow_2(local_energy(lattice,conf,loc)/2);
  }
  free(loc);
  return(energy2/conf.elements);
}
