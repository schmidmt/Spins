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
#include <blocking.h>

int
main (int argc, char **argv)
{
  settings conf,conf1,conf2;
  static int specified_file_flag = 0;
  static int verbose_flag = 0;
  char conf_file[PATH_MAX + NAME_MAX];
  char outputfile_name[PATH_MAX + NAME_MAX];
  int c; /* Output status for getopt*/
  int useclud,random_spins;
  int rng_seed;
  int i,j,k;
  int step,steps_of_beta;
  double beta, beta_start,beta_end,beta_step;
  FILE * outputfp;
  datapoint data;

  const gsl_rng_type * RngType;

  gsl_rng_env_setup();
  RngType = gsl_rng_default;
  conf.rng = gsl_rng_alloc (RngType);

  lattice_site * lattice   = NULL;
  lattice_site * lattice1  = NULL;
  lattice_site * lattice2  = NULL;

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

  /*Copy old settings to blocked settings*/
  memcpy(&conf1,&conf,sizeof(settings));
  conf1.sidelength = conf.sidelength/2;
  conf1.elements   = conf.elements/4;
  memcpy(&conf2,&conf1,sizeof(settings));
  conf2.sidelength = conf1.sidelength/2;
  conf2.elements   = conf1.elements/4;

  printf("#Allocating\n");
  lattice  = allocate_lattice(conf);
  lattice1 = allocate_lattice(conf1);
  lattice2 = allocate_lattice(conf2);

  randomize_spins(lattice,conf);
  
  /**********************
   * Running Simulation *
   **********************/
  beta          = beta_start;
  steps_of_beta = (int) ceil((beta_end-beta_start)/beta_step);
  step          = 0;
  //                   1     2      3     4     5     6    7
  fprintf(outputfp,"#%-11s %-13s %-13s %-13s %-13s %-13s %-13s\n","Beta","<m>","m err","<E>","E err","Specific Heat","Magnetic Susc");
  printf("#Starting simulation\n\n");

  double * e_block, * m_block, * e_block_avg, * m_block_avg, \
         * e_block_error, * m_block_error, * c_block , * chi_block;
  gsl_vector * mag_vector;

  e_block       = (double *) malloc(conf.block_size*sizeof(double));
  m_block       = (double *) malloc(conf.block_size*sizeof(double));
  e_block_avg   = (double *) malloc(conf.blocks*sizeof(double));
  m_block_avg   = (double *) malloc(conf.blocks*sizeof(double));
  c_block       = (double *) malloc(conf.blocks*sizeof(double));
  chi_block     = (double *) malloc(conf.blocks*sizeof(double));
  e_block_error = (double *) malloc(conf.blocks*sizeof(double));
  m_block_error = (double *) malloc(conf.blocks*sizeof(double));
  mag_vector = gsl_vector_alloc(conf.spindims);

  for(k = 0 ; k < steps_of_beta ; k++)
  {
    loadBar(step,steps_of_beta,50,80);
    fflush(stdout);

    /*****************************************************************/

    //Settle first
    for(i = 0 ; i < conf.max_settle ; i++)
    {
     clusterupdate(lattice,conf,beta);
    }

    //Get averages and stdev for messurements
    for(i = 0 ; i < conf.blocks ; i++)
    {
      for(j = 0 ; j < conf.block_size ; j++)
      {
        clusterupdate(lattice,conf,beta);
        block_4_majority(conf,lattice,conf1,lattice1);
        //block_4_majority(conf1,lattice1,conf2,lattice2);
        e_block[j] = total_energy(lattice1,conf1);
        m_block[j] = magnetization(lattice1,conf1,mag_vector);
      }
      e_block_avg[i]   = gsl_stats_mean(e_block,1,conf.block_size);
      e_block_error[i] = gsl_stats_sd(e_block,1,conf.block_size);
      m_block_avg[i]   = gsl_stats_mean(m_block,1,conf.block_size);
      m_block_error[i] = gsl_stats_sd(m_block,1,conf.block_size);
      c_block[i]       = gsl_pow_2(beta)*gsl_pow_2(e_block_error[i]);
      chi_block[i]     = beta*gsl_pow_2(m_block_error[i]);
    }
    data.beta      = beta;
    data.erg       = gsl_stats_mean(e_block_avg,1,conf.blocks);
    data.erg_error = gsl_stats_sd(e_block_avg,1,conf.blocks);
    data.mag       = gsl_stats_mean(m_block_avg,1,conf.blocks);
    data.mag_error = gsl_stats_sd(m_block_avg,1,conf.blocks);
    data.c         = gsl_stats_mean(c_block,1,conf.blocks);
    data.c_error   = gsl_stats_sd(c_block,1,conf.blocks);
    data.chi       = gsl_stats_mean(chi_block,1,conf.blocks);
    data.chi_error = gsl_stats_sd(chi_block,1,conf.blocks);
    /*****************************************************************/

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
