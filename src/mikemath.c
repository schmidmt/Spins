#include <../src/mikemath.h>

int intpow(int x, int y)
{
  /*
  int ans = 1;
  if(y == 0) return 1;
  if(y == 1) return x;
  while( y > 0 )
  { 
    ans *= x;
    y--;
  }
  */
  return((int)pow(x,y));
}

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
