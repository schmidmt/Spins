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
#include <mikemath.h>
#include <lattice.h>
#include <physics.h>

int
main (int argc, char **argv)
{
  static int verbose_flag;
  static int specified_file_flag = 0;
  char conf_file[PATH_MAX + NAME_MAX];
  int c; /* Output status for getopt*/
  int useclud,random_spins;
  int spacedims, spindims, sidelength;
  int screenwidth,rng_seed;
  int i,j,k;
  int printiter,max_settle;
  double beta, beta_start,beta_end,beta_step;
  double * mag_record;
  double mag_mean, mag_stdev, mag_mean_last, mag_stdev_last;
  int mag_record_size;
  /*screenwidth = atoi(getenv("COLUMNS"));*/

  const gsl_rng_type * RngType;
  gsl_rng * rng;

  gsl_rng_env_setup();
  RngType = gsl_rng_default;
  rng = gsl_rng_alloc (RngType);

  gsl_vector **lattice = NULL;

  int * location = (int *) malloc(spacedims*sizeof(int));
  int * neigh    = (int *) malloc(spacedims*sizeof(int));


  /* Physical Quantities */
  double mag;


  /* Libconf Stuff */
  config_t cfg;
  config_setting_t *setting;
  const char *str;

  config_init(&cfg);

  /* Processing arguments */
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


  /* Reading Config File */
  if(!config_read_file(&cfg,conf_file))
  {
    fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),config_error_line(&cfg), config_error_text(&cfg));
    config_destroy(&cfg);
    return(EXIT_FAILURE);
  }

  for(i = 0 ; i < 10 ; i++)
    printf("#");
  printf("\n");

  if(config_lookup_string(&cfg, "name", &str))
    printf("#%-27s = '%s'\n","Run Name",str);
  else
  {
    fprintf(stderr,"#No 'name' setting found in configuration file.\n");
    return(EXIT_FAILURE);
  }
  if(config_lookup_int(&cfg,"spindims",&spindims))
    printf("#%-27s = %d\n","Number of spin dimensions",spindims);
  else
  {
    printf("#No Specified spindims, defaulting to 2\n");
    spindims = 2;
  }
  if(config_lookup_int(&cfg,"spacedims",&spacedims))
    printf("#%-27s = %d\n","Number of space dimensions",spacedims);
  else
  {
    printf("#No Specified spacedims, defaulting to 2\n");
    spacedims = 2;
  }
  if(config_lookup_int(&cfg,"sidelength",&sidelength))
    printf("#%-27s = %d\n","Sidelenth of lattice",sidelength);
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
    printf("#%-27s = %s\n","Random Spins","True");
  else
    printf("#%-27s = %s\n","Random Spins","False");
    

  if(!config_lookup_bool(&cfg,"useclud",&useclud))
  {
    printf("#No Specified update algorithm, defaulting to metropolis\n");
    useclud = 0;
  }
  if(useclud != 0)
    printf("#%-27s = %s\n","Update Algorith","Cluster");
  else
    printf("#%-27s = %s\n","Update Algorith","Metropolis");

  if(config_lookup_float(&cfg,"beta_start",&beta_start))
    printf("#%-27s = %e\n","Beta Start",beta_start);
  else
  {
    printf("#No Specified beta_start, defaulting to 10\n");
    beta_start = 10;
  }

  if(config_lookup_float(&cfg,"beta_end",&beta_end))
    printf("#%-27s = %e\n","Beta End",beta_end);
  else
  {
    printf("#No Specified beta_end, defaulting to 10\n");
    beta_end = 10;
  }
  if(config_lookup_float(&cfg,"beta_step",&beta_step))
    printf("#%-27s = %e\n","Beta Step",beta_end);
  else
  {
    printf("#No Specified beta_step, defaulting to 10\n");
    beta_step= 0.001;
  }


  if(config_lookup_int(&cfg,"printiter",&printiter))
    printf("#%-27s = %d\n","Printiter",printiter);
  else
  {
    printf("#No Specified printiter, defaulting to 10\n");
    printiter = 10;
  }

  if(config_lookup_int(&cfg,"max_settle",&max_settle))
    printf("#%-27s = %d\n","Max Settle Iterations",max_settle);
  else
  {
    printf("#No Specified max_settle, defaulting to 1000\n");
    max_settle= 10000;
  }
  if(config_lookup_int(&cfg,"rng_seed",&rng_seed))
  {
    printf("#%-27s = %d\n","Random Seed",rng_seed);
    gsl_rng_set(rng,rng_seed);
  }


  /* mag_record_size is set to the size of the lattice.
   * We'll see if that is enough to see a settling.
   */
  /*mag_record_size = div(max_settle,20).quot;*/
  mag_record_size = intpow(sidelength,spacedims);
  mag_record = (double *) calloc(mag_record_size,sizeof(double));

  for(i = 0 ; i < 10 ; i++)
    printf("#");
  printf("\n");


  /* ALLOCATION OF LATTICE AND VECTORS */
  
  lattice = allocate_lattice(sidelength,spacedims,spindims);

  if(verbose_flag)
    fprintf(stderr,"#Allocated %d points on the lattice\n",intpow(sidelength,spacedims));
  
  
  if(random_spins == 1)
    randomize_spins(lattice,sidelength,spacedims,spindims,rng);
  else
    set_homogenious_spins(lattice,sidelength,spacedims,spindims);

  if(verbose_flag) 
    print_lattice (lattice,sidelength,spacedims,spindims);

  
  /* Running Simulation */
  beta = beta_start;
  while(beta <= beta_end)
  {
    mag_mean_last  = DBL_MAX;
    mag_stdev_last = DBL_MAX;
    for(i = 0 ; i < mag_record_size ; i++)
      mag_record[i] = 0.0;

    for(i = 0 ; i < max_settle ; i++)
    {
      metropolis_update(lattice,sidelength,spacedims,spindims,beta,rng);
      mag = magnetization(lattice,sidelength,spacedims,spindims);
      enqueue(mag_record,mag_record_size,mag);
      if(i > mag_record_size && i%(mag_record_size) == 0)
      {
        mag_mean  = gsl_stats_mean(mag_record,1,mag_record_size);
        mag_stdev = gsl_stats_sd(mag_record,1,mag_record_size);
        /* Consider settled if stdev and mean is within 10% */
        if( fabs((mag_mean-mag_mean_last)/mag_mean_last) < 0.1 && fabs((mag_stdev-mag_stdev_last)/mag_stdev_last) < 0.1)
        {
          /*fprintf(stderr,"Stopped after i/lattice size = %d",i/intpow(sidelength,spacedims));*/
          break;
        }
        mag_mean_last  = mag_mean;
        mag_stdev_last = mag_stdev;
      }
      /*
      printf("%05d %e %e %e\n",i,mag,mag_mean+2*mag_stdev,mag_mean-2*mag_stdev);
      */
    }
    printf("%e %e\n",beta,mag_mean);
    beta += beta_step;
  }

  /* CLEANUP */
  free_lattice(lattice,sidelength,spacedims);
  free(location);
  location = NULL;
  free(neigh);
  neigh = NULL;
  free(mag_record);
  mag_record = NULL;

  /* Free GSL Random Num Generator*/
  gsl_rng_free (rng);


  return EXIT_SUCCESS;
}
