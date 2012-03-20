#ifndef COMMON_H
#define COMMON_H 1

#include <math.h>
#include <float.h>
#include <gsl/gsl_rng.h>

#define intpow(x,y) (int)pow(x,y)

typedef struct
{
  int sidelength;
  int spacedims;
  int spindims;
  int useclud;
  char * outputfile_name;
  char * conf_file;
  int verbose_flag;
  int max_settle;
  gsl_rng * rng;
  int elements;
  int block_size;
  int blocks;
} settings;

typedef struct
{
  double beta;
  double erg;
  double erg_error;
  double mag;
  double mag_error;
} datapoint;




void
enqueue(double * array, int size, double new);

inline void
loadBar(int x, int n, int r, int w);


#endif /* !COMMON_H */
