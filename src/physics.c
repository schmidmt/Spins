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
magnetization(lattice_site * lattice, settings conf, gsl_vector * mag_vector)
{
  int i;
  gsl_vector_set_zero(mag_vector);
  double result = 0;
  for(i = 0 ; i < conf.elements ; i++)
  {
    gsl_vector_add(mag_vector,lattice[i].spin);
  }
  gsl_blas_ddot(mag_vector,mag_vector,&result);

  /* I have no idea why this is happening! */
  return(result/intpow(conf.elements,2));
  //return(result);
}


double
local_energy (lattice_site * lattice, settings conf, int loc_id)
{
  double energy = 0,result;
  int i,secondary;

  for(i = 0 ; i < 2*(conf.spacedims) ; i++)
  {
    secondary = lattice[loc_id].neighbors[i];
    gsl_blas_ddot(lattice[loc_id].spin,lattice[secondary].spin,&result);
    energy -= result;
  }
  //return(energy/(intpow(sidelength,spacedims)));
  return(energy);
}


double
new_local_energy (lattice_site * lattice, settings conf, int loc_id, gsl_vector * new)
{
  double energy = 0,result;
  int i,secondary;
  
  for(i = 0 ; i < 2*(conf.spacedims) ; i++)
  {
    secondary = lattice[loc_id].neighbors[i];
    gsl_blas_ddot(new,lattice[secondary].spin,&result);
    energy -= result;
  }
  return(energy);
}

double
total_energy(lattice_site * lattice, settings conf)
{
  double energy = 0,result;
  int i,j,secondary;

  for(j = 0 ; j < conf.elements ; j++)
  {
    for(i = 0 ; i < 2*(conf.spacedims) ; i += 2)
    {
      secondary = lattice[j].neighbors[i];
      gsl_blas_ddot(lattice[j].spin,lattice[secondary].spin,&result);
      energy -= result;
    }
  }
  
  return(energy/conf.elements);
}
