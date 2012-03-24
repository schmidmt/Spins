/******************************************************************************
  Copyright 2012 Michael Schmidt (mts@colorado.edu)

  This file is part of spins.
 *
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
 *******************************************************************************/


#include <stdio.h>
#include <string.h>
#include <lattice.h>
#include <gsl/gsl_rng.h>
#include <common.h>

int
main (int argc, char * argv[])
{
  int i,j;
  int size = argc - 2;
  int *data  = (int *)  malloc(sizeof(int)*size);
  int sidelength, spacedims;
  int * loc;
  int * neigh;
  int num;
  gsl_rng * rng;
  const gsl_rng_type * RngType;
  gsl_rng_env_setup();
  RngType = gsl_rng_default;
  rng = gsl_rng_alloc (RngType);

  /* Read in data */
  if(size != 0)
  {
    for (i = 0 ; i< size ; i++)
    {
      data[i] = atoi(argv[i+2]);
    }
  }
  switch(atoi(argv[1]))
  {
    case 0: /* Neighbor */
      sidelength = data[0];
      spacedims  = data[1];
      loc   = (int *) malloc(sizeof(int)*spacedims);
      neigh = (int *) malloc(sizeof(int)*spacedims);
      for(i = 0; i < spacedims ; i++)
      {
        loc[i] = data[2+i];
      }
      for(i = 0; i < 2*spacedims ; i++)
      {
        neighbor(sidelength,spacedims, loc, neigh, i);
        for(j = 0; j < spacedims ; j++)
        {
          printf("%d ",neigh[j]);
        }
        
         if( i != 2*spacedims-1)
          printf(", ");
      }
      printf("\n");
      break;
    case 1: /* NUM_TO_LOCATION */
      sidelength = data[0];
      spacedims  = data[1];
      num        = data[2];
      loc   = (int *) malloc(sizeof(int)*spacedims);
      num_to_location(sidelength,spacedims,num,loc);
      for(i = 0 ; i < spacedims ; i++)
      {
        printf("%d",loc[i]);
        if(i != spacedims-1)
          printf(" ");
      }
      printf("\n");
      break;
    case 2: /* location_to_num */
      sidelength = data[0];
      spacedims  = data[1];
      loc   = (int *) malloc(sizeof(int)*spacedims);
      for(i = 0; i < spacedims ; i++)
      {
        loc[i] = data[2+i];
      }
      num =  location_to_num(sidelength,spacedims,loc);
      printf("%d\n",num);
      break;
    case 3: /* num to location -> location to num */
      sidelength = data[0];
      spacedims  = data[1];
      loc   = (int *) malloc(sizeof(int)*spacedims);
      for(i = 0 ; i < 100 ; i++)
      {
        num = gsl_rng_uniform_int(rng,intpow(sidelength,spacedims)-1);
        num_to_location(sidelength,spacedims,num,loc);
        if( num != location_to_num(sidelength,spacedims,loc) )
        {
          printf("FAILURE\n");
          exit(EXIT_FAILURE);
        }
      }
      printf("SUCCESS\n");
      break;
    default:
      printf("No arguments!\n");
      exit(EXIT_FAILURE);
  }

}
