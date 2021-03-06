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


/* clusterupdate.c: Tools for updating a lattice using the Wolff Algorithm
 */

#include <stdio.h>
#include <stdlib.h>
#include <common.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <lattice.h>
#include <math.h>
#include <clusterupdate.h>
#include <physics.h>

/********************************************************************************
 * clusterupdatebatch: Runs clusterupdate multiple times and gets physics as well
 * as error estimates.
 *******************************************************************************/
int
clusterupdatebatch(lattice_site * lattice, settings conf, double beta, datapoint * data )
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
      e_block[j] = total_energy(lattice,conf);
      m_block[j] = magnetization(lattice,conf,mag_vector);
    }
    e_block_avg[i]   = gsl_stats_mean(e_block,1,conf.block_size);
    e_block_error[i] = gsl_stats_sd(e_block,1,conf.block_size);
    m_block_avg[i]   = gsl_stats_mean(m_block,1,conf.block_size);
    m_block_error[i] = gsl_stats_sd(m_block,1,conf.block_size);
    c_block[i]       = gsl_pow_2(beta)*gsl_pow_2(e_block_error[i]);
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
/*******************************************************************************
 * clusterupdate: Runs the Wolff algorithm on lattice.
 ******************************************************************************/
int
clusterupdate(lattice_site * lattice, settings conf, double beta)
{
  int i,j,cluster_size, update_start;
  gsl_vector * base, * delta;
  int * update_list;
  int ** bonds;
  double scale;

  base  = gsl_vector_alloc(conf.spindims); /* Direction to switch over */
  delta = gsl_vector_alloc(conf.spindims); /* How much to change a vector */
  update_list = (int *) calloc(conf.elements,sizeof(int));   /* List of Lattice points to switch */
  bonds       = (int **) malloc(conf.elements*sizeof(int*)); /* Bonds Matrix to make sure bonds aren't rechecked. */

  for(i = 0 ; i < conf.elements ; i++)
    bonds[i]  = (int *) calloc(conf.elements,sizeof(int));

  /* Generate a random unit vector to set as base */
  unit_vec(base,conf.rng);

  //Pick a random point on the lattice then pass off to gencluster
  update_start = gsl_rng_uniform_int(conf.rng,conf.elements);
  update_list[update_start] = 1;
  
  //Pass off to gencluster
  cluster_size = gencluster(lattice,conf,update_start,update_list,bonds,base,beta);
  cluster_size++; //So it will include the first element
  
  //Flip the entire cluster
  for(j = 0 ; j < conf.elements ; j++)
  {
    if(update_list[j] == 1)
    {
      gsl_blas_ddot(lattice[j].spin,base,&scale);
      gsl_vector_memcpy(delta,base);
      gsl_vector_scale(delta,-2.0*scale);
      gsl_vector_add(lattice[j].spin,delta);
    }
  }

  gsl_vector_free(base);
  gsl_vector_free(delta);
  free(update_list);
  for(i = 0 ; i < conf.elements ; i++)
    free(bonds[i]);
  free(bonds);
  return(cluster_size);
}

/*******************************************************************************
 * gencluster: recursivily calls itself and returns a set of lattice points.
 ******************************************************************************/
int
gencluster(lattice_site * lattice, settings conf, int loc_id , int * update_list, int ** bonds, gsl_vector * base, double beta)
{
  int i,update_count = 0;
  double exp_factor,s1n,s2n;
  int nid;

  for(i = 0 ; i < 2*conf.spacedims ; i++)
  {
    nid = lattice[loc_id].neighbors[i]; /* Neighbor ID */
    //Continue on if the point has already been added
    //  or if bond has already been checked.
    if(update_list[nid] == 1 || bonds[loc_id][nid] != 0)
      continue;
    gsl_blas_ddot(lattice[loc_id].spin,base,&s1n);
    gsl_blas_ddot(lattice[nid].spin,base,&s2n);
    exp_factor = -2.0*s1n*s2n*beta;
  
    if(exp_factor > -10 && gsl_rng_uniform(conf.rng) < 1-gsl_sf_exp(exp_factor))
    {
      update_list[nid] = 1; /* Set to be flipped */
      update_count += 1;    
      bonds[loc_id][nid] = 1; /* Prevent link from being checked again */
      bonds[nid][loc_id] = 1;
      update_count += gencluster(lattice,conf,nid,update_list,bonds,base,beta);
    }
    else
    {
      bonds[loc_id][nid] = -1; /* If not added keep this link from being made again */
      bonds[nid][loc_id] = -1;
    }
  }
  return(update_count);
}
