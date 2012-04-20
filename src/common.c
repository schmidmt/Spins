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


#include <common.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>

void
enqueue(double * array, int size, double new)
{
  /* enqueue shifts every value left and sets
   * the last element to new.
   */
  int i;
  for(i = 0 ; i < size-1 ; i++)
  {
    array[i] = array[i+1];
  }
  array[size-1] = new;
}


// Process has done x out of n rounds,
// and we want a bar of width w and resolution r.
inline void
loadBar(int x, int n, int r, int w)
{
  if(r > n || r == 0)
  {
    r = n;
  }
  // Only update r times.
  if ( x % (n/r) != 0 )
  {
    return;
  }
  printf("\33[1A\33[2K");

  // Calculuate the ratio of complete-to-incomplete.
  float ratio = x/(float)n;
  int   c     = ratio * w;

  // Show the percentage complete.
  printf("%3d%% [", (int)(ratio*100) );

  // Show the load bar.
  for (x=0; x<c; x++)
    printf("=");

  for (x=c; x<w; x++)
    printf(" ");

  printf("]\n");
}

void
print_data(FILE * fh, datapoint data)
{
  fprintf(fh,"%+e %+e %+e %+e %+e %+e %+e %+e %+e\n",data.beta,data.mag,data.mag_error,data.erg,data.erg_error,data.c,data.c_error,data.chi,data.chi_error);
}


inline void
unit_vec(gsl_vector * vect , gsl_rng * rng )
{
  double sqrsum = 0, random_num;
  int j;

  for(j = 0 ; j < vect->size ; j++)
  {
    random_num = 2*(gsl_rng_uniform(rng)-0.5);
    gsl_vector_set(vect,j,random_num);
    sqrsum += gsl_pow_2(random_num);
  }
  gsl_vector_scale(vect,1.0/sqrt(sqrsum));
}

void
print_vec(gsl_vector * vect, char * name)
{
  int i;
  printf("%s = [ ",name);
  i = 0;
  if( (int)vect->size > 1 )
  {
    for(i = 0 ; i < ((int)vect->size-1) ; i++)
    {
      printf("%+e, ",gsl_vector_get(vect,i));
    }
    i++;
  }
  printf("%+e ]\n",gsl_vector_get(vect,i));
}
