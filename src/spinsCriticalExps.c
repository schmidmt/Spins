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
#include <metropolis.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <common.h>
#include <lattice.h>
#include <physics.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_math.h>

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
  double beta,beta_start,beta_end,beta_delta,beta_step;
  int beta_data_points;
  datapoint data;

  double * tau_log, * c_log, *c_error_log, * chi_log, 
         * chi_error_log, *m_error_log, *m_log;
  
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
  printf("#Starting simulation\n\n");
  beta_start  = CRITB+0.0001;
  beta_delta  = 0.0001;
  beta_end    = CRITB+0.100;

  beta_data_points = (int)floor((beta_end-beta_start)/beta_delta);

  tau_log       = (double *) malloc(beta_data_points*sizeof(double));
  c_log         = (double *) malloc(beta_data_points*sizeof(double));
  c_error_log   = (double *) malloc(beta_data_points*sizeof(double));
  m_log         = (double *) malloc(beta_data_points*sizeof(double));
  m_error_log   = (double *) malloc(beta_data_points*sizeof(double));
  chi_log       = (double *) malloc(beta_data_points*sizeof(double));
  chi_error_log = (double *) malloc(beta_data_points*sizeof(double));

  beta = beta_start;

  for(i = 0 ; i < beta_data_points ; i++)
  {
    loadBar(i,beta_data_points,50,80);
    //clusterupdatebatch(lattice,conf,beta,&data);
    mupdatebatch(lattice,conf,beta,&data);
    tau_log[i]       = gsl_sf_log(fabs(1-CRITB/beta));
    c_log[i]         = gsl_sf_log(data.c);
    c_error_log[i]   = 1.0/gsl_pow_2(gsl_sf_log(data.c_error));
    m_log[i]         = gsl_sf_log(data.mag);
    m_error_log[i]   = 1.0/gsl_pow_2(gsl_sf_log(fabs(data.mag_error)));
    chi_log[i]       = gsl_sf_log(data.chi);
    chi_error_log[i] = 1.0/gsl_pow_2(gsl_sf_log(data.chi_error));

    fprintf(outputfp,"%g %g %g\n",tau_log[i],chi_log[i],c_log[i]);
    beta += beta_delta;
  }

  double c0,c1,cov00,cov01,cov11,chisq;

  gsl_fit_wlinear(tau_log,1,c_error_log,1,c_log,1,beta_data_points,
                  &c0, &c1, &cov00, &cov01, &cov11, &chisq);
  printf("# alpha = %g with chisq = %g\n",-c1,chisq/intpow(beta_data_points,2));
  gsl_fit_wlinear(tau_log,1,m_error_log,1,m_log,1,beta_data_points,
                  &c0, &c1, &cov00, &cov01, &cov11, &chisq);
  printf("# beta = %g with chisq = %g\n",-c1,chisq/intpow(beta_data_points,2));
  gsl_fit_wlinear(tau_log,1,chi_error_log,1,chi_log,1,beta_data_points,
                  &c0, &c1, &cov00, &cov01, &cov11, &chisq);
  printf("# gamma = %g with chisq = %g\n",-c1,chisq/intpow(beta_data_points,2));

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
