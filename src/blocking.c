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
block_5_majority(settings confin,lattice_site * lin , settings confout, \
                 lattice_site * lout)
{
  int i,j,k;
  if(confin.spindims != 1)
  {
    fprintf(stderr,"Cannot run majority rule on more than 1 spin dimension!\n");
    return(-1);
  }

  

  return(0);
}




