#ifndef MIKEMATH_H
#define MIKEMATH_H 1

#include <math.h>
#include <float.h>

/*
int
intpow(int x, int y);
*/

#define intpow(x,y) (int)pow(x,y)


void
enqueue(double * array, int size, double new);

inline void
loadBar(int x, int n, int r, int w);


#endif /* !MIKEMATH_H */
