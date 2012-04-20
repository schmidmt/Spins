#ifndef CLUSTERUPDATE_H
#define CLUSTERUPDATE_H 1

#include <gsl/gsl_vector.h>
#include <common.h>

int gencluster(lattice_site * lattice, settings conf, int loc_id , int * update_list, int ** bonds, gsl_vector * base, double beta);
int clusterupdate(lattice_site * lattice, settings conf, double beta);
int clusterupdatebatch(lattice_site * lattice, settings conf, double beta, datapoint * data );


#endif /* !CLUSTERUPDATE_H */
