#include <../src/mikemath.h>

int intpow(int x, int y)
{
  int ans = 1;
  if(y == 0) return 1;
  if(y == 1) return x;
  while( y > 0 )
  { 
    ans *= x;
    y--;
  }
  return(ans);
}
