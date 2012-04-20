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
  int i;
  int step,steps_of_beta;
  double beta, beta_start,beta_end,beta_step;
  FILE * outputfp;
  datapoint data;

  const gsl_rng_type * RngType;

  gsl_rng_env_setup();
  RngType = gsl_rng_default;
  conf.rng = gsl_rng_alloc (RngType);

  lattice_site * lattice   = NULL;

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
  lattice = allocate_lattice(conf);

  randomize_spins(lattice,conf);
  
  /**********************
   * Running Simulation *
   **********************/
  beta          = beta_start;
  steps_of_beta = (int) ceil((beta_end-beta_start)/beta_step);
  step          = 0;
  //                   1     2      3     4     5     6    7
  fprintf(outputfp,"#%-13s %-13s %-13s %-13s %-13s %-13s %-13s\n","Beta","<m>","m err","<E>","E err","Specific Heat","Magnetic Susc");

  printf("#Starting simulation\n\n");
  while(step <= steps_of_beta)
  {
    if(!conf.verbose_flag)
    {
      loadBar(step,steps_of_beta,50,80);
    }
    fflush(stdout);

    mupdatebatch(lattice,conf,beta,&data);
    print_data(outputfp,data);
    beta += beta_step;
    step++;
  }

  /***********
   * CLEANUP *
   ***********/
  printf("#Cleaning Up\n");
  free_lattice(lattice,conf);
  fclose(outputfp);
  gsl_rng_free (conf.rng);
  printf("#Done\n");
  return EXIT_SUCCESS;
}
