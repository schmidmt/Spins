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
  static int verbose_flag;
  static int specified_file_flag = 0;
  char conf_file[PATH_MAX + NAME_MAX];
  char outputfile_name[PATH_MAX + NAME_MAX];
  int c; /* Output status for getopt*/
  int useclud,random_spins;
  int spacedims, spindims, sidelength;
  int rng_seed;
  int i,j;
  int max_settle;
  int step,steps_of_beta;
  double beta, beta_start,beta_end,beta_step;
  /*screenwidth = atoi(getenv("COLUMNS"));*/
  FILE * outputfp;

  const gsl_rng_type * RngType;
  gsl_rng * rng;

  gsl_rng_env_setup();
  RngType = gsl_rng_default;
  rng = gsl_rng_alloc (RngType);

  gsl_vector **lattice = NULL;
  gsl_vector * magnet  = NULL ;

  int * location = (int *) malloc(spacedims*sizeof(int));
  int * neigh    = (int *) malloc(spacedims*sizeof(int));

  /* Physical Quantities */
  double mag;
  double energy;
  double ** data;


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

  if (verbose_flag)
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

  if(config_lookup_int(&cfg,"spindims",&spindims))
  {
    printf("#%-27s = %d\n","Number of spin dimensions",spindims);
    fprintf(outputfp,"#%-27s = %d\n","Number of spin dimensions",spindims);
  }
  else
  {
    printf("#No Specified spindims, defaulting to 2\n");
    spindims = 2;
  }

  if(config_lookup_int(&cfg,"spacedims",&spacedims))
  {
    printf("#%-27s = %d\n","Number of space dimensions",spacedims);
    fprintf(outputfp,"#%-27s = %d\n","Number of space dimensions",spacedims);
  }
  else
  {
    printf("#No Specified spacedims, defaulting to 2\n");
    spacedims = 2;
  }

  if(config_lookup_int(&cfg,"sidelength",&sidelength))
  {
    printf("#%-27s = %d\n","Sidelenth of lattice",sidelength);
    fprintf(outputfp,"#%-27s = %d\n","Sidelenth of lattice",sidelength);
  }
  else
  {
    printf("#No Specified sidelength, defaulting to 10\n");
    sidelength = 10;
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
    beta_start = 0;
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

  if(config_lookup_int(&cfg,"max_settle",&max_settle))
  {
    printf("#%-27s = %d\n","Max Settle Iterations",max_settle);
    fprintf(outputfp,"#%-27s = %d\n","Max Settle Iterations",max_settle);
  }
  else
  {
    printf("#No Specified max_settle, defaulting to 1000\n");
    max_settle= 10000;
  }

  if(config_lookup_int(&cfg,"rng_seed",&rng_seed))
  {
    printf("#%-27s = %d\n","Random Seed",rng_seed);
    fprintf(outputfp,"#%-27s = %d\n","Random Seed",rng_seed);
    gsl_rng_set(rng,rng_seed);
  }

  for(i = 0 ; i < 10 ; i++)
  {
    printf("#");
    fprintf(outputfp,"#");
  }
  printf("\n");
  fprintf(outputfp,"\n");
  fprintf(outputfp,"#%-11s %-12s %-12s\n","Beta","<m>","<E>");


  /*************************************
   * ALLOCATION OF LATTICE AND VECTORS *
   *************************************/
  printf("#Allocating\n");
  
  lattice = allocate_lattice(sidelength,spacedims,spindims);

  if(verbose_flag)
    fprintf(stderr,"#Allocated %d points on the lattice\n",intpow(sidelength,spacedims));
  
  
  if(random_spins == 1)
  {
    if(verbose_flag)
      fprintf(stderr,"Randomizing Spins\n");
    randomize_spins(lattice,sidelength,spacedims,spindims,rng);
  }
  else
  {
    if(verbose_flag)
      fprintf(stderr,"Setting Homogenious Spins\n");
    set_homogenious_spins(lattice,sidelength,spacedims,spindims);
    //set_checkerboard_spins(lattice,sidelength,spacedims,spindims);
  }

  /* if(verbose_flag) 
    print_lattice (lattice,sidelength,spacedims,spindims); */

  magnet = gsl_vector_alloc(spindims);

  
  /**********************
   * Running Simulation *
   **********************/
  beta          = beta_start;
  mag           = magnetization(lattice,sidelength,spacedims,spindims,magnet);
  energy        = total_energy(lattice,sidelength,spacedims,spindims);
  steps_of_beta = (int) ceil((beta_end-beta_start)/beta_step);
  step          = 0;

  data = (double **) malloc(sizeof(double **)*(steps_of_beta+1)); 
  for(i = 0 ; i <= steps_of_beta ; i ++)
  {
    data[i] = (double *) calloc(3,sizeof(double));
  }
 
  printf("#Starting simulation\n\n");
  while(step <= steps_of_beta)
  {
    data[step][0] = beta; 
    loadBar(step,steps_of_beta,50,80);
    fflush(stdout);

    /* This specifies how man steps the
     * particular beta will be evaluated at 
     * if it never settles down.
     */
    for(i = 0 ; i < max_settle ; i++)
    {
      metropolis_update(lattice,sidelength,spacedims,spindims,beta,rng,magnet,&energy);
    }
    gsl_blas_ddot(magnet,magnet,&mag);
    energy = total_energy(lattice,sidelength,spacedims,spindims);
    mag = magnetization(lattice,sidelength,spacedims,spindims,magnet);
    fprintf(outputfp,"%e %e %e\n",beta,mag/intpow(sidelength,spacedims),energy/intpow(sidelength,spacedims));
    data[step][1] = mag;
    data[step][2] = energy;
    beta += beta_step;
    step++;
  }
  printf("\n");

  /***********
   * CLEANUP *
   ***********/
  printf("#Cleaning Up\n");
  free_lattice(lattice,sidelength,spacedims);
  free(location);
  location = NULL;
  free(neigh);
  neigh = NULL;
  fclose(outputfp);

  /* Free GSL Random Num Generator*/
  gsl_rng_free (rng);
  
  printf("#Done\n");

  return EXIT_SUCCESS;
}
