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
  static int specified_file_flag = 0;
  static int verbose_flag = 0;
  char conf_file[PATH_MAX + NAME_MAX];
  char outputfile_name[PATH_MAX + NAME_MAX];
  int c; /* Output status for getopt*/
  int useclud,random_spins;
  int rng_seed;
  int i,j;
  int step,steps_of_beta;
  double beta, beta_start,beta_end,beta_step;
  FILE * outputfp;


  const gsl_rng_type * RngType;

  gsl_rng_env_setup();
  RngType = gsl_rng_default;
  conf.rng = gsl_rng_alloc (RngType);

  gsl_vector ** lattice    = NULL;
  gsl_vector * mag_vector  = NULL ;


  /* Physical Quantities */
  double mag,mag_error,mag2,mag2_error;
  double energy,energy_error,energy2,energy2_error;
  double magsus,specheat;

  /* Libconf Stuff */
  config_t cfg;
  //config_setting_t *setting;
  const char *str;

  config_init(&cfg);


  /*******************************************************************
   * This file includes all of the argv and conf file parsing stuff. *
   *******************************************************************/
  #include "conf_file_parser.c"

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
  mag           = magnetization(lattice,conf,mag_vector);
  energy        = total_energy(lattice,conf);
  steps_of_beta = (int) ceil((beta_end-beta_start)/beta_step);
  step          = 0;
  //                   1     2      3     4     5     6    7
  fprintf(outputfp,"#%-11s %-13s %-13s %-13s %-13s %-13s %-13s\n","Beta","<m>","m err","<E>","E err","Specific Heat","Magnetic Susc");

  //mupdate(lattice,conf,beta,mag_vector,&mag,&mag_error,&energy,&energy_error);
  //print_lattice(lattice,conf);
  printf("#Starting simulation\n\n");
  while(step <= steps_of_beta)
  {
    if(!conf.verbose_flag)
    {
      loadBar(step,steps_of_beta,50,80);
    }
    fflush(stdout);

    mupdate(lattice,conf,beta,mag_vector,&mag,&mag_error,&mag2,&mag2_error,&energy,&energy_error,&energy2,&energy2_error);
    gsl_blas_ddot(mag_vector,mag_vector,&mag);
    energy = total_energy(lattice,conf);
    mag = magnetization(lattice,conf,mag_vector);
    specheat = (energy2 - gsl_pow_2(energy/conf.elements));
    magsus   = mag2 - gsl_pow_2(mag/conf.elements);
    fprintf(outputfp,"%e %+e %+e %+e %+e %+e %+e\n",beta,fabs(mag)/conf.elements,mag_error/conf.elements,energy/conf.elements,energy_error/conf.elements,specheat/conf.elements,magsus/conf.elements);
    beta += beta_step;
    step++;
  }
  //printf("\n");
  //
  //print_lattice(lattice,conf);

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
