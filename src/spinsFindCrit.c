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
  conf.sidelength = 8;
  conf.max_settle  = 100;
  conf.block_size = 100;
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
  double var_l, var_c, var_r;
  beta_start      = 0.4;
  beta_step       = 0.0001;
  double dbeta;
  beta = beta_start;
  while(err > 0.0001)
  {
    //Find first derivative
    mupdatebatch(lattice,conf,beta-(beta_step/2.0),&data);
    var_l = data.c;
    mupdatebatch(lattice,conf,beta+(beta_step/2.0),&data);
    var_r = data.c;
    firstd = (var_r-var_l)/beta_step; 

    //Find second derivative
    mupdatebatch(lattice,conf,beta-beta_step,&data);
    var_l = data.c;
    mupdatebatch(lattice,conf,beta,&data);
    var_c = data.c;
    mupdatebatch(lattice,conf,beta+beta_step,&data);
    var_r = data.c;
    secondd = (var_l+var_r-2*var_c)/pow(beta_step,2); 
    
    dbeta = firstd/secondd;
    err = fabs(dbeta);
    beta = beta - dbeta;
    printf("Beta Est = %g  Err = %g\n",beta,err);
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
