#ifndef COMMON_H
#define COMMON_H 1

#include <math.h>
#include <float.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>

#define intpow(x,y) (int)pow(x,y)
#define CRITT 1.0/0.44
#define CRITB 0.44

typedef struct
{
  int sidelength;
  int spacedims;
  int spindims;
  int useclud;
  char * outputfile_name;
  char * conf_file;
  int output2file;
  int verbose_flag;
  int max_settle;
  gsl_rng * rng;
  int elements;
  int block_size;
  int blocks;
} settings;

typedef struct
{
  gsl_vector * spin;
  int * neighbors;
} lattice_site;

typedef struct
{
  double beta;
  double erg;
  double erg_error;
  double mag;
  double mag_error;
  double c;
  double c_error;
  double chi;
  double chi_error;
} datapoint;




void
enqueue(double * array, int size, double new);

inline void
loadBar(int x, int n, int r, int w);

void
print_data( FILE * fh, datapoint data);

inline void
unit_vec(gsl_vector * vect, gsl_rng * rng);

inline void
print_vec(gsl_vector * vect, char * name);

#endif /* !COMMON_H */
