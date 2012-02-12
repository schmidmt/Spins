#include <../src/mikemath.h>
#include <stdio.h>
#include <stdlib.h>


int
main (void)
{
  int answer = 0;
  // 2^0 Case
  answer = intpow(2,0);
  printf("intpow(2,0) = %d\n",answer);
  if(answer != 1)
    return(EXIT_FAILURE);

  // 2^1 Case
  answer = intpow(2,1);
  printf("intpow(2,1) = %d\n",answer);
  if(answer != 2)
    return(EXIT_FAILURE);

  // 2^10 Case
  answer = intpow(2,10);
  printf("intpow(2,10) = %d\n",answer);
  if(answer != 1024)
    return(EXIT_FAILURE);

  // 3^5 Case
  answer = intpow(3,5);
  printf("intpow(3,5) = %d\n",answer);
  if(answer != 243)
    return(EXIT_FAILURE);

  return(EXIT_SUCCESS);
}
