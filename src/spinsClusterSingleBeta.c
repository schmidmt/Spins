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
  static int specified_file_flag = 0;
  static int verbose_flag = 0;
  char conf_file[PATH_MAX + NAME_MAX];
  char outputfile_name[PATH_MAX + NAME_MAX];
  int c; /* Output status for getopt*/
  int useclud,random_spins;
  int rng_seed;
  int i;
  int step,steps_of_beta;
  double beta, beta_start,beta_end,beta_step;
  FILE * outputfp;


  const gsl_rng_type * RngType;

  gsl_rng_env_setup();
  RngType = gsl_rng_default;
  conf.rng = gsl_rng_alloc (RngType);

  lattice_site * lattice   = NULL;
  gsl_vector * mag_vector  = NULL ;

  /* Libconf Stuff */
  config_t cfg;
  //config_setting_t *setting;
  const char *str;

  config_init(&cfg);

  /*******************************************************************
  * This file includes all of the argv and conf file parsing stuff. *
  *******************************************************************/
  #include <conf_file_parser.c>

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
  
  
  if(random_spins == 1)
  {
    if(conf.verbose_flag)
      fprintf(stderr,"Randomizing Spins\n");
    randomize_spins(lattice,conf);
  }
  else
  {
    if(conf.verbose_flag)
      fprintf(stderr,"Setting Homogenious Spins\n");
    set_homogenious_spins(lattice,conf);
    //set_checkerboard_spins(lattice,sidelength,conf.spacedims,conf.spindims);
  }

  /* if(conf.verbose_flag) 
    print_lattice (lattice,sidelength,conf.spacedims,conf.spindims); */

  mag_vector = gsl_vector_calloc(conf.spindims);
  
  /**********************
   * Running Simulation *
   **********************/
  beta          = beta_start;
  steps_of_beta = (int) ceil((beta_end-beta_start)/beta_step);
  step          = 0;
  //                   1     2      3     4     5     6    7
  fprintf(outputfp,"#%-11s %-13s %-13s %-13s %-13s %-13s %-13s\n","Beta","<m>","m err","<E>","E err","Specific Heat","Magnetic Susc");
  
  /*
  //print_lattice(lattice,conf);
  printf("%d %e %e\n",-1,total_energy(lattice,conf)/conf.elements,magnetization(lattice,conf,mag_vector)/conf.elements);
  for(i = 0 ; i < 400 ; i++)
  {
    clusterupdate(lattice,conf,beta);
    //print_lattice(lattice,conf);
    printf("%d %e %e\n",i,total_energy(lattice,conf),magnetization(lattice,conf,mag_vector));
  }
  */
  datapoint data;
  data.beta = 0;
  clusterupdatebatch(lattice,conf,beta,&data);
  printf("%e %e %e %e %e\n",data.beta,data.erg,data.erg_error,data.mag,data.mag_error);

  /***********
   * CLEANUP *
   ***********/
  printf("#Cleaning Up\n");
  free_lattice(lattice,conf);
  free(location);
  location = NULL;
  free(neigh);
  neigh = NULL;
  fclose(outputfp);

  /* Free GSL Random Num Generator*/
  gsl_rng_free (conf.rng);
  

  printf("#Done\n");

  return EXIT_SUCCESS;
}
