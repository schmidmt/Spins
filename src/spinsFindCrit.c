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


#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <getopt.h>
#include <libconfig.h>
#include <limits.h>
#include <string.h>
#include <metropolis.h>
#include <clusterupdate.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <common.h>
#include <lattice.h>
#include <physics.h>

int
main (int argc, char **argv)
{
  settings conf;
  double beta_step, beta_start, beta;
  datapoint data;
  double c0, cl1, cl2, cr1, cr2;

  const gsl_rng_type * RngType;

  gsl_rng_env_setup();
  RngType = gsl_rng_default;
  conf.rng = gsl_rng_alloc (RngType);

  lattice_site * lattice   = NULL;
  gsl_vector * mag_vector  = NULL ;


  /************
   * SETTINGS *
   ************/
  conf.spindims   = 2;
  conf.spacedims  = 2;
  conf.sidelength = 32;
  conf.max_settle = 1000;
  conf.block_size = 1000;
  conf.blocks     = 20;
  conf.verbose_flag = 0;


  /*************************************
   * ALLOCATION OF LATTICE AND VECTORS *
   *************************************/
  //Store Element number for future use
  conf.elements = intpow(conf.sidelength,conf.spacedims);

  printf("#Allocating\n");


  int * location = (int *) malloc(conf.spacedims*sizeof(int));
  int * neigh    = (int *) malloc(conf.spacedims*sizeof(int));
  
  lattice = allocate_lattice(conf);

  if(conf.verbose_flag)
    fprintf(stderr,"#Allocated %d points on the lattice\n",conf.elements);
  
  
  randomize_spins(lattice,conf);

  /* if(conf.verbose_flag) 
    print_lattice (lattice,sidelength,conf.spacedims,conf.spindims); */

  mag_vector = gsl_vector_calloc(conf.spindims);
  
  /**********************
   * Running Simulation *
   **********************/
  double err = 100;
  double firstd, secondd;
  beta_start      = 0.45;
  beta_step       = 0.001;
  double dbeta;
  beta = beta_start;
  while(err > 0.01)
  {
    //Calculate values of specific heat
    //c(beta)
    clusterupdatebatch(lattice,conf,beta,&data);
    c0 = data.c;
    //Left and Right one step
    clusterupdatebatch(lattice,conf,beta-beta_step,&data);
    cl1 = data.c;
    clusterupdatebatch(lattice,conf,beta+beta_step,&data);
    cr1 = data.c;
    //Left and Right one step
    clusterupdatebatch(lattice,conf,beta-2*beta_step,&data);
    cl2 = data.c;
    clusterupdatebatch(lattice,conf,beta+2*beta_step,&data);
    cr2 = data.c;

    //Find first derivative
    firstd  = ( ((cl2-cr2)/12) +(2*(-cl1+cr1)/3) )/beta_step;
    //Find second derivative
    secondd = ( ((-cl2-cr2)/12) +(4*(cl1+cr1)/3) - (5*c0/2) )/pow(beta_step,2);
    
    dbeta = firstd/secondd;
    err = fabs(firstd);
    beta = beta - dbeta;
    printf("Beta Est = %g  delta= %g firstd= %g secondd= %g\n",beta,dbeta,firstd,secondd);
  }

  /***********
   * CLEANUP *
   ***********/
  printf("#Cleaning Up\n");
  free_lattice(lattice,conf);
  free(location);
  location = NULL;
  free(neigh);
  neigh = NULL;

  /* Free GSL Random Num Generator*/
  gsl_rng_free (conf.rng);
  

  printf("#Done\n");

  return EXIT_SUCCESS;
}
