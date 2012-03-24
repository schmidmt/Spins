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


#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include <lattice.h>
#include <common.h>
#include <stdlib.h>
#include <string.h>

gsl_vector **
allocate_lattice (settings conf)
{
  int i;
  gsl_vector ** lattice = (gsl_vector **) malloc(conf.elements*sizeof(gsl_vector *));
  for ( i = 0 ; i < conf.elements ; i++ )
  {
    lattice[i] =  (gsl_vector *) gsl_vector_alloc (conf.spindims);
  }
  return(lattice);
}


int
free_lattice ( gsl_vector ** lattice , settings conf )
{
  int i;
  for ( i = 0 ; i < conf.elements ; i++ )
  {
    gsl_vector_free(lattice[i]);
  }
  free(lattice);
  lattice = NULL;
  return 0;
}

void
print_lattice (gsl_vector ** lattice, settings conf )
{
  int i,j;
  int * location = (int *) malloc(conf.spacedims*sizeof(int));
  for(i = 0 ; i < conf.elements ; i++ )
  {
    location = num_to_location(conf,i,location);
    printf("lattice[%d]:  (",i);
    for(j = 0 ; j < conf.spacedims ; j++)
    {
      printf(" %d",location[j]);
    }
    printf(" ) -> ( ");
    for(j = 0 ; j < conf.spindims ; j++)
      printf("%e ",gsl_vector_get (lattice[i], j));
    printf(")\n");
  }
  free(location);
}

void
set_homogenious_spins(gsl_vector ** lattice, settings conf )
{
  int i;
  for(i = 0 ; i < conf.elements ; i++)
  {
    gsl_vector_set_basis(lattice[i],0);
  }
}

void
set_checkerboard_spins(gsl_vector ** lattice, settings conf)
{
  int i;
  for(i = 0 ; i < conf.elements ; i++)
  {
    gsl_vector_set_basis(lattice[i],0);
    if(i % 2 == 1)
    {
      gsl_vector_scale(lattice[i],-1.0);
    }
  }
}

void
randomize_spins(gsl_vector ** lattice, settings conf )
{
  int i,j;
  double sqrsum,random_num;
  
  for(i = 0 ; i < conf.elements ; i++)
  {
    sqrsum = 0;
    for(j = 0 ; j < conf.spindims ; j++)
    {
      random_num = 2*(gsl_rng_uniform(conf.rng)-0.5);
      gsl_vector_set(lattice[i],j,random_num);
      sqrsum += gsl_pow_2(random_num);
    }
    gsl_vector_scale(lattice[i],1.0/sqrt(sqrsum));
  }
}

int 
neighbor(settings conf, int * loc, int * neigh , int num)
{
  int dir = 0, ind = 0;
  if((num & 1) == 0)
    dir = -1;
  else
    dir = 1;
  
  ind = div(num,2).quot;
  if(ind >= conf.spacedims)
  {
    fprintf(stderr,"Trying to find index out of range...\n");
    exit(EXIT_FAILURE);
  }

  memcpy(neigh,loc,conf.spacedims*sizeof(int));

  neigh[ind] += dir;
  if(neigh[ind] < 0)
    neigh[ind] += conf.sidelength;
  if(neigh[ind] >= conf.sidelength)
    neigh[ind] += -(conf.sidelength);

  return(0);
}

int *
num_to_location(settings conf, int num, int * location)
{
  int i,power;
  //int * location = (int *) malloc(spacedims*sizeof(int));
  div_t asdf;
  for(i = 0 ; i < conf.spacedims ; i++)
  {
    power = intpow(conf.sidelength,conf.spacedims-1-i);
    asdf = div(num,power);
    location[i] = asdf.quot;
    num = asdf.rem;
  }

  return(location);
}

int
location_to_num ( settings conf, int * location )
{
  int num = 0;
  int i;
  for(i = 0 ; i < conf.spacedims ; i++)
  {
    num += location[i]*intpow(conf.sidelength,conf.spacedims-1-i);
  }
  return(num);
}
