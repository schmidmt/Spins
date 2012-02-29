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
  double mag,mag_error;
  double energy,energy_error;


  /* Libconf Stuff */
  config_t cfg;
  //config_setting_t *setting;
  const char *str;

  config_init(&cfg);

  /************************
   * Processing arguments *
   ************************/
  while (1)
    {
      static struct option long_options[] = {
        {"verbose", no_argument, &verbose_flag, 1},
        {"brief", no_argument, &verbose_flag, 0},
        {"file", required_argument, 0, 'f'},
        {0, 0, 0, 0}
      };
      /* getopt_long stores the option index here. */
      int option_index = 0;

      c = getopt_long (argc, argv, "f:", long_options, &option_index);

      if (c == -1)
        break;

      switch (c)
      {
        case 0:
          if (long_options[option_index].flag != 0)
            break;
          printf ("option %s", long_options[option_index].name);
          if (optarg)
            printf (" with arg %s", optarg);
          printf ("\n");
          break;

        case 'f':
          printf ("#Using Config file: '%s'\n", optarg);
          strcpy(conf_file,optarg);
          specified_file_flag = 1;
          break;

        case '?':
          /* getopt_long already printed an error message. */
          break;

        default:
          abort ();
      }
    }

  conf.verbose_flag = verbose_flag;

  if (conf.verbose_flag)
    puts ("#verbose flag is set");

  /* Print any remaining command line arguments (not options). */
  if (optind < argc)
    {
      printf ("non-option ARGV-elements: ");
      while (optind < argc)
	{
	  printf ("%s ", argv[optind++]);
	}
      putchar ('\n');
    }


  if(!specified_file_flag)
  {
    fprintf(stderr,"#No file specified, please specify one with the -f flag\n");
    exit(EXIT_FAILURE);
  }


  /***********************
   * Reading Config File *
   ***********************/
  if(!config_read_file(&cfg,conf_file))
  {
    fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),config_error_line(&cfg), config_error_text(&cfg));
    config_destroy(&cfg);
    return(EXIT_FAILURE);
  }

  if(config_lookup_string(&cfg, "name", &str))
  {
    strcpy(outputfile_name,str);
    strcat(outputfile_name,".dat");
    if((outputfp = fopen(outputfile_name,"w+")) == NULL )
    {
      fprintf(stderr,"Cannot open file %s\n",outputfile_name);
      exit(EXIT_FAILURE);
    }
    for(i = 0 ; i < 10 ; i++)
    {
      printf("#");
      fprintf(outputfp,"#");
    }
    printf("\n");
    fprintf(outputfp,"\n");
    printf("#%-27s = '%s'\n","Run Name",str);
    fprintf(outputfp,"#%-27s = '%s'\n","Run Name",str);
  }
  else
  {
    fprintf(stderr,"#No 'name' setting found in configuration file.\n");
    return(EXIT_FAILURE);
  }

  if(config_lookup_int(&cfg,"spindims",&conf.spindims))
  {
    printf("#%-27s = %d\n","Number of spin dimensions",conf.spindims);
    fprintf(outputfp,"#%-27s = %d\n","Number of spin dimensions",conf.spindims);
  }
  else
  {
    printf("#No Specified spindims, defaulting to 2\n");
    conf.spindims = 2;
  }

  if(config_lookup_int(&cfg,"spacedims",&conf.spacedims))
  {
    printf("#%-27s = %d\n","Number of space dimensions",conf.spacedims);
    fprintf(outputfp,"#%-27s = %d\n","Number of space dimensions",conf.spacedims);
  }
  else
  {
    printf("#No Specified spacedims, defaulting to 2\n");
    conf.spacedims = 2;
  }

  if(config_lookup_int(&cfg,"sidelength",&conf.sidelength))
  {
    printf("#%-27s = %d\n","Sidelenth of lattice",conf.sidelength);
    fprintf(outputfp,"#%-27s = %d\n","Sidelenth of lattice",conf.sidelength);
  }
  else
  {
    printf("#No Specified sidelength, defaulting to 10\n");
    conf.sidelength = 10;
  }

  if(!config_lookup_bool(&cfg,"randomize_spins",&random_spins))
  {
    printf("#Randomize_spins not specified, defaulting to random.\n");
    random_spins = 0;
  }

  if(random_spins == 1)
  {
    printf("#%-27s = %s\n","Random Spins","True");
    fprintf(outputfp,"#%-27s = %s\n","Random Spins","True");
  }
  else
  {
    printf("#%-27s = %s\n","Random Spins","False");
    fprintf(outputfp,"#%-27s = %s\n","Random Spins","False");
  }

  if(!config_lookup_bool(&cfg,"useclud",&useclud))
  {
    printf("#No Specified update algorithm, defaulting to metropolis\n");
    useclud = 0;
  }

  if(useclud != 0)
  {
    printf("#%-27s = %s\n","Update Algorith","Cluster");
    fprintf(outputfp,"#%-27s = %s\n","Update Algorith","Cluster");
  }
  else
  {
    printf("#%-27s = %s\n","Update Algorith","Metropolis");
    fprintf(outputfp,"#%-27s = %s\n","Update Algorith","Metropolis");
  }

  if(config_lookup_float(&cfg,"beta_start",&beta_start))
  {
    printf("#%-27s = %e\n","Beta Start",beta_start);
    fprintf(outputfp,"#%-27s = %e\n","Beta Start",beta_start);
  }
  else
  {
    printf("#No Specified beta_start, defaulting to 0\n");
    beta_start = 1;
  }

  if(config_lookup_float(&cfg,"beta_end",&beta_end))
  {
    printf("#%-27s = %e\n","Beta End",beta_end);
    fprintf(outputfp,"#%-27s = %e\n","Beta End",beta_end);
  }
  else
  {
    printf("#No Specified beta_end, defaulting to 10\n");
    beta_end = 10;
  }

  if(config_lookup_float(&cfg,"beta_step",&beta_step))
  {
    printf("#%-27s = %e\n","Beta Step",beta_end);
    fprintf(outputfp,"#%-27s = %e\n","Beta Step",beta_end);
  }
  else
  {
    printf("#No Specified beta_step, defaulting to 0.1\n");
    beta_step= 0.1;
  }

  if(config_lookup_int(&cfg,"max_settle",&conf.max_settle))
  {
    printf("#%-27s = %d\n","Max Settle Iterations",conf.max_settle);
    fprintf(outputfp,"#%-27s = %d\n","Max Settle Iterations",conf.max_settle);
  }
  else
  {
    printf("#No Specified max_settle, defaulting to 1000\n");
    conf.max_settle= 10000;
  }

  if(config_lookup_int(&cfg,"rng_seed",&rng_seed))
  {
    printf("#%-27s = %d\n","Random Seed",rng_seed);
    fprintf(outputfp,"#%-27s = %d\n","Random Seed",rng_seed);
    gsl_rng_set(conf.rng,rng_seed);
  }
  if(config_lookup_int(&cfg,"block_size",&conf.block_size))
  {
    printf("#%-27s = %d\n","Blocking Size",conf.block_size);
    fprintf(outputfp,"#%-27s = %d\n","Blocking Size",conf.block_size);
  }
  else
  {
    printf("#No Specified block_size, defaulting to 1000\n");
    conf.block_size = 1000;
  }
  if(config_lookup_int(&cfg,"blocks",&conf.blocks))
  {
    printf("#%-27s = %d\n","Blocks",conf.blocks);
    fprintf(outputfp,"#%-27s = %d\n","Blocks",conf.blocks);
  }
  else
  {
    printf("#No Specified blocks, defaulting to 10\n");
    conf.blocks = 10;
  }

  for(i = 0 ; i < 10 ; i++)
  {
    printf("#");
    fprintf(outputfp,"#");
  }
  printf("\n");
  fprintf(outputfp,"\n");


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
  fprintf(outputfp,"#%-11s %-13s %-13s %-13s %-13s\n","Beta","<m>","m err","<E>","E err");

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

    mupdate(lattice,conf,beta,mag_vector,&mag,&mag_error,&energy,&energy_error);
    //gsl_blas_ddot(mag_vector,mag_vector,&mag);
    //energy = total_energy(lattice,conf);
    //mag = magnetization(lattice,conf,mag_vector);
    fprintf(outputfp,"%e %+e %+e %+e %+e\n",beta,fabs(mag)/conf.elements,mag_error/conf.elements,energy/conf.elements,energy_error/conf.elements);
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
