#ifndef CLUSTERUPDATE_H
#define CLUSTERUPDATE_H 1

#include <gsl/gsl_vector.h>
#include <common.h>

int gencluster(gsl_vector ** lattice, settings conf, int * loc , int * update_list, gsl_vector * base, double beta);
int clusterupdate(gsl_vector ** lattice, settings conf, double beta);
int clusterupdatebatch(gsl_vector ** lattice, settings conf, double beta, datapoint * data );


#endif /* !CLUSTERUPDATE_H */
