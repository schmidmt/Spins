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
