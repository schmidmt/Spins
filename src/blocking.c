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


int
block_4_majority(settings confin,lattice_site * lin , settings confout, \
                 lattice_site * lout)
{
  int i,j;
  int * newloc, * oldloc;
  int oid, nid;
  double newspin;

  if(confin.spindims != 1)
  {
    fprintf(stderr,"Cannot run majority rule on more than 1 spin dimension!\n");
    return(-1);
  }
  if((confin.sidelength & 1) != 0)
  {
    fprintf(stderr,"Cannot run on lattices of odd sidelength!\n");
    return(-1);
  }

  if(confout.sidelength != confin.sidelength/2)
  {
    fprintf(stderr,"confout isn't set correctly.\n");
    return(-1);
  }

  oldloc = malloc(confin.spacedims*sizeof(double));
  newloc = malloc(confout.spacedims*sizeof(double));

  for(i = 0 ; i < confout.sidelength; i++)
  {
    newloc[0] = i;
    for(j = 0 ; j < confout.sidelength; j++)
    {
      newspin = 0;
      newloc[1] = j;
      nid = location_to_num(confout,newloc);
      
      /* Lower Left */
      oldloc[0] = 2*i;
      oldloc[1] = 2*j;
      oid = location_to_num(confin,oldloc);
      newspin += gsl_vector_get(lin[oid].spin,0);

      /* Upper Left */
      oldloc[0] = 2*i;
      oldloc[1] = 2*j+1;
      oid = location_to_num(confin,oldloc);
      newspin += gsl_vector_get(lin[oid].spin,0);

      /* Lower Right */
      oldloc[0] = 2*i+1;
      oldloc[1] = 2*j;
      oid = location_to_num(confin,oldloc);
      newspin += gsl_vector_get(lin[oid].spin,0);

      /* Upper Right */
      oldloc[0] = 2*i+1;
      oldloc[1] = 2*j+1;
      oid = location_to_num(confin,oldloc);
      newspin += gsl_vector_get(lin[oid].spin,0);

      if(newspin > 0)
        gsl_vector_set(lout[nid].spin,0,1);
      else if(newspin < 0)
        gsl_vector_set(lout[nid].spin,0,-1);
      else if( gsl_rng_uniform(confout.rng) > 0.5)
        gsl_vector_set(lout[nid].spin,0,-1);
      else
        gsl_vector_set(lout[nid].spin,0,1);
    }
  }
  
  free(newloc);
  free(oldloc);
  return(0);
}
