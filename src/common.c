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
