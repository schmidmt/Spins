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


#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>
#include <math.h>
#include <lattice.h>
#include <common.h>
#include <physics.h>

int
mupdate_step(gsl_vector ** lattice, settings conf, double beta, gsl_vector * mag_vector, double * energy )
{
  if(conf.verbose_flag)
    printf("metropolis_update: ");
  int i;
  int * loc = (int *) malloc(conf.spacedims*sizeof(int));
  double de,dep,sqrsum,random_num,exp_factor;
  gsl_vector * newspin = gsl_vector_alloc(conf.spindims);
  int change_flag = 0;

  /* Create a substitute random spin direction
   * which we will compare to our present spin
   */
  sqrsum = 0;
  for(i = 0 ; i < conf.spindims ; i++)
  {
    random_num = 2*(gsl_rng_uniform(conf.rng)-0.5);
    gsl_vector_set(newspin,i,random_num);
    sqrsum += gsl_pow_2(random_num);
  }
  gsl_vector_scale(newspin,1.0/sqrt(sqrsum));

  if(conf.verbose_flag)
  {
    printf("newspin = ( ");
    for(i = 0 ; i < conf.spindims; i++)
    {
      printf("%+1.1e ",gsl_vector_get(newspin,i));
    }
    printf(") ");
    printf("loc = ( ");
  }
  /* Choose a random location on the lattice */
  for(i = 0 ; i < conf.spacedims ; i++)
  {
    loc[i] = gsl_rng_uniform_int(conf.rng,conf.sidelength);
    if(conf.verbose_flag)
      printf("%d ",loc[i]);
  }
  if(conf.verbose_flag)
  {
    printf(") ");

    printf("oldspin = ( ");
    for(i = 0 ; i < conf.spindims; i++)
    {
      printf("%+1.1e ",gsl_vector_get(lattice[location_to_num(conf,loc)],i));
    }
    printf(") ");
  }
  
  de  = local_energy(lattice,conf,loc);
  dep = new_local_energy(lattice,conf,loc,newspin);

  if(conf.verbose_flag)
    printf("de = %+1.1e dep = %+1.1e ",de,dep);

  exp_factor = -beta*(dep-de);

  /* Calculate the probibility of acceptance of the new vector, pacc.
   * If the expontential facor is less than -10, it's basically zero,
   * this code prevents it from throwing an underflow exception.
   */

  if(dep < de)  
    change_flag = 1;
  else if(exp_factor > -10 && gsl_rng_uniform(conf.rng) < gsl_sf_exp(exp_factor))
  {
    change_flag = 1;
  }

  /* Generate a uniform random number between (0,1],
   * if it's < pacc, switch the vector.
   */
  if(change_flag)
  {
    //gsl_vector_sub(mag_vector,lattice[location_to_num(sidelength,spacedims,loc)]);
    //gsl_vector_add(mag_vector,newspin);
    //*energy += (dep-de);
    gsl_vector_memcpy(lattice[location_to_num(conf,loc)],newspin);
    if(conf.verbose_flag)
      printf("Y\n");
  }
  else
  {
    if(conf.verbose_flag)
      printf("N\n");
  }
  free(loc);
  gsl_vector_free(newspin);
  newspin = NULL;
  return(0);
}

/* mupdate: Thermalizes lattice and calculates errors.
 */
int
mupdate(gsl_vector ** lattice, settings conf, double beta, gsl_vector * mag_vector, double * mag, double * mag_error, double * mag2, double * mag2_error, double * energy , double * energy_error, double * energy2, double* energy2_error)
{
  int i,j;
  double * energy_log = NULL, * mag_log = NULL;
  double * energy_mean = NULL, * mag_mean = NULL;
  double * energy_err_log = NULL, * mag_err_log = NULL;
  double * energy2_log = NULL, * mag2_log = NULL;
  double * energy2_mean = NULL, * mag2_mean = NULL;
  double * energy2_err_log = NULL, * mag2_err_log = NULL;

  energy_log      = (double *) malloc(sizeof(double *)*(conf.block_size));
  energy_mean     = (double *) malloc(sizeof(double *)*(conf.blocks));
  energy_err_log  = (double *) malloc(sizeof(double *)*(conf.blocks));
  mag_log         = (double *) malloc(sizeof(double *)*(conf.block_size));
  mag_mean        = (double *) malloc(sizeof(double *)*(conf.blocks));
  mag_err_log     = (double *) malloc(sizeof(double *)*(conf.blocks));
  energy2_log     = (double *) malloc(sizeof(double *)*(conf.block_size));
  energy2_mean    = (double *) malloc(sizeof(double *)*(conf.blocks));
  energy2_err_log = (double *) malloc(sizeof(double *)*(conf.blocks));
  mag2_log        = (double *) malloc(sizeof(double *)*(conf.block_size));
  mag2_mean       = (double *) malloc(sizeof(double *)*(conf.blocks));
  mag2_err_log    = (double *) malloc(sizeof(double *)*(conf.blocks));
  //Thermalize 
  /*
  FILE * therm_out;
  therm_out = fopen("thermout.dat","w+");
  */
  for(i = 0 ; i < conf.max_settle ; i++)
  {
    mupdate_step(lattice, conf, beta, mag_vector, energy);
    //fprintf(therm_out,"%d %e %e\n",i,total_energy(lattice,conf)/conf.elements,magnetization(lattice,conf,mag_vector)/conf.elements);
  }
  //fclose(therm_out);
  //return(0);

  /**************************************
   * Run more to get averages and error *
   **************************************/
  for(j = 0 ; j < conf.blocks ; j++) //Number of blocks
  {
    for(i = 0 ; i < conf.block_size ; i++) //Blocksize
    {
      mupdate_step(lattice, conf, beta, mag_vector, energy );
      energy_log[i] = total_energy(lattice,conf);
      mag_log[i]    = magnetization(lattice,conf,mag_vector);
      energy2_log[i] = total_energy2(lattice,conf);
      mag2_log[i]    = magnetization2(lattice,conf,mag_vector);
    }
    energy_mean[j]    = gsl_stats_mean(energy_log,1,conf.block_size);
    energy_err_log[j] = gsl_stats_sd(energy_log,1,conf.block_size);
    mag_mean[j]       = gsl_stats_mean(mag_log,1,conf.block_size);
    mag_err_log[j]    = gsl_stats_sd(mag_log,1,conf.block_size);
    energy2_mean[j]    = gsl_stats_mean(energy2_log,1,conf.block_size);
    energy2_err_log[j] = gsl_stats_sd(energy2_log,1,conf.block_size);
    mag2_mean[j]       = gsl_stats_mean(mag2_log,1,conf.block_size);
    mag2_err_log[j]    = gsl_stats_sd(mag2_log,1,conf.block_size);
  }
  *energy       = gsl_stats_mean(energy_mean,1,conf.blocks);
  *energy_error = gsl_stats_mean(energy_err_log,1,conf.blocks);
  *mag          = gsl_stats_mean(mag_mean,1,conf.blocks);
  *mag_error    = gsl_stats_mean(mag_err_log,1,conf.blocks);
  *energy2       = gsl_stats_mean(energy_mean,1,conf.blocks);
  *energy2_error = gsl_stats_mean(energy_err_log,1,conf.blocks);
  *mag2          = gsl_stats_mean(mag_mean,1,conf.blocks);
  *mag2_error    = gsl_stats_mean(mag_err_log,1,conf.blocks);

  free(energy_log);
  energy_log = NULL;
  free(energy_mean);
  energy_mean = NULL;
  free(energy_err_log);
  energy_err_log = NULL;
  free(mag_log);
  mag_log = NULL;
  free(mag_mean);
  mag_mean = NULL;
  free(mag_err_log);
  mag_err_log = NULL;
  free(energy2_log);
  energy2_log = NULL;
  free(energy2_mean);
  energy2_mean = NULL;
  free(energy2_err_log);
  energy2_err_log = NULL;
  free(mag2_log);
  mag2_log = NULL;
  free(mag2_mean);
  mag2_mean = NULL;
  free(mag2_err_log);
  mag2_err_log = NULL;
  //Compute Errors
  return(0);
}
