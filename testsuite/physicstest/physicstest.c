#include <stdio.h>
#include <string.h>
#include <lattice.h>
#include <physics.h>
#include <gsl/gsl_rng.h>
#include <common.h>

int
main (int argc, char * argv[])
{
  int i,j;
  int size = argc - 2;
  int *data  = (int *)  malloc(sizeof(int)*size);
  int sidelength, spacedims, spindims;
  int * loc;
  int num;
  gsl_rng * rng;
  const gsl_rng_type * RngType;
  gsl_rng_env_setup();
  RngType = gsl_rng_default;
  rng = gsl_rng_alloc (RngType);
  gsl_vector ** lattice;
  gsl_vector * magnet;
  double mag,energy;

  /* Read in data */
  if(size != 0)
  {
    for (i = 0 ; i< size ; i++)
    {
      data[i] = atoi(argv[i+2]);
    }
  }
  switch(atoi(argv[1]))
  {
    case 0: /* Magnetization */
      /**********************************************
       * Outputs magnetization of a uniform lattice *
       **********************************************/
      sidelength = data[0];
      spacedims  = data[1];
      spindims   = data[2];
      magnet = gsl_vector_alloc(spindims);
      lattice = allocate_lattice(sidelength,spacedims,spindims);
      set_homogenious_spins(lattice,sidelength,spacedims,spindims);
      mag = magnetization(lattice,sidelength,spacedims,spindims,magnet);
      free_lattice(lattice,sidelength,spacedims);
      printf("%2.1f\n",mag);
      break;
    case 1: /* Local Energy */
      /**********************************************
       * Outputs energy of a uniform lattice point *
       **********************************************/
      sidelength = data[0];
      spacedims  = data[1];
      spindims   = data[2];
      loc        = (int *) malloc(sizeof(int)*spacedims);
      for(i = 0 ; i < spacedims ; i++)
        loc[i] = data[i+3];
      lattice = allocate_lattice(sidelength,spacedims,spindims);
      set_homogenious_spins(lattice,sidelength,spacedims,spindims);
      energy = 0;
      energy = local_energy(lattice, sidelength, spacedims, spindims, loc);
      free_lattice(lattice,sidelength,spacedims);
      printf("%1.3e\n",energy);
      break;
    case 2: /* Total Energy */
      /********************************************
       * Outputs energy of a checkerboard lattice *
       ********************************************/
      sidelength = data[0];
      spacedims  = data[1];
      spindims   = data[2];
      lattice = allocate_lattice(sidelength,spacedims,spindims);
      set_checkerboard_spins(lattice,sidelength,spacedims,spindims);
      energy = total_energy(lattice, sidelength, spacedims, spindims );
      printf("%1.3e\n",energy);
      break;
    default:
      printf("No arguments!\n");
      exit(EXIT_FAILURE);
  }

}
