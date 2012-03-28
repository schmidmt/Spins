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
mupdate_step(lattice_site * lattice, settings conf, double beta)
{
  int i;
  int loc_id;
  double de,dep,sqrsum,random_num,exp_factor;
  gsl_vector * newspin = gsl_vector_alloc(conf.spindims);

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

  /* Choose a random location on the lattice */
  loc_id = gsl_rng_uniform_int(conf.rng,conf.elements);
  
  de  = local_energy(lattice,conf,loc_id);
  dep = new_local_energy(lattice,conf,loc_id,newspin);

  exp_factor = -beta*(dep-de);

  /* Calculate the probibility of acceptance of the new vector, pacc.
   * If the expontential facor is less than -10, it's basically zero,
   * this code prevents it from throwing an underflow exception.
   */

  if(dep < de || (exp_factor > -10 && gsl_rng_uniform(conf.rng) < gsl_sf_exp(exp_factor)))
  {
    gsl_vector_memcpy(lattice[loc_id].spin,newspin);
  }

  gsl_vector_free(newspin);
  newspin = NULL;
  return(0);
}

/* mupdate: Thermalizes lattice and calculates errors.
 */
int
mupdatebatch(lattice_site * lattice, settings conf, double beta, datapoint * data)
{
  int i,j;
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

  //Settle First
  for(i = 0 ; i < conf.max_settle ; i++)
  {
    mupdate_step(lattice, conf, beta);
  }

  /**************************************
   * Run more to get averages and error *
   **************************************/
  for(i = 0 ; i < conf.blocks ; i++) //Number of blocks
  {
    for(j = 0 ; j < conf.block_size ; j++) //Blocksize
    {
      mupdate_step(lattice, conf, beta);
      e_block[j] = total_energy(lattice,conf);
      m_block[j] = magnetization(lattice,conf,mag_vector);
    }
    e_block_avg[i]   = gsl_stats_mean(e_block,1,conf.block_size);
    e_block_error[i] = gsl_stats_sd(e_block,1,conf.block_size);
    m_block_avg[i]   = gsl_stats_mean(m_block,1,conf.block_size);
    m_block_error[i] = gsl_stats_sd(m_block,1,conf.block_size);
    c_block[i]       = beta*gsl_pow_2(e_block_error[i]);
    chi_block[i]     = beta*gsl_pow_2(m_block_error[i]);
  }
  (*data).beta      = beta;
  (*data).erg       = gsl_stats_mean(e_block_avg,1,conf.blocks);
  (*data).erg_error = gsl_stats_sd(e_block_avg,1,conf.blocks);
  (*data).mag       = gsl_stats_mean(m_block_avg,1,conf.blocks);
  (*data).mag_error = gsl_stats_sd(m_block_avg,1,conf.blocks);
  (*data).c         = gsl_stats_mean(c_block,1,conf.blocks);
  (*data).c_error   = gsl_stats_sd(c_block,1,conf.blocks);
  (*data).chi       = gsl_stats_mean(chi_block,1,conf.blocks);
  (*data).chi_error = gsl_stats_sd(chi_block,1,conf.blocks);

  free(e_block);
  free(m_block);
  free(e_block_avg);
  free(m_block_avg);
  free(e_block_error);
  free(m_block_error);
  gsl_vector_free(mag_vector);
  return(0);
}
