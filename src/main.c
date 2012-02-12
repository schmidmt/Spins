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
#include <mikemath.h>
#include <lattice.h>

int
main (int argc, char **argv)
{
  static int verbose_flag;
  static int specified_file_flag = 0;
  char conf_file[PATH_MAX + NAME_MAX];
  int c; /* Output status for getopt*/
  int useclud;
  int spacedims, spindims, sidelength;
  int screenwidth;
  int i,j,k;
  /*screenwidth = atoi(getenv("COLUMNS"));*/

  gsl_vector **lattice = NULL;

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
	  printf ("Using Config file: '%s'\n", optarg);
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
    puts ("verbose flag is set");

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
    fprintf(stderr,"No file specified, please specify one with the -f flag\n");
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
    printf("-");
  printf("\n");

  if(config_lookup_string(&cfg, "name", &str))
    printf("%-27s = '%s'\n","Run Name",str);
  else
  {
    fprintf(stderr,"No 'name' setting found in configuration file.\n");
    return(EXIT_FAILURE);
  }
  if(config_lookup_int(&cfg,"spindims",&spindims))
    printf("%-27s = %d\n","Number of spin dimensions",spindims);
  else
  {
    printf("No Specified spindims, defaulting to 2\n");
    spindims = 2;
  }
  if(config_lookup_int(&cfg,"spacedims",&spacedims))
    printf("%-27s = %d\n","Number of space dimensions",spacedims);
  else
  {
    printf("No Specified spacedims, defaulting to 2\n");
    spacedims = 2;
  }
  if(config_lookup_int(&cfg,"sidelength",&sidelength))
    printf("%-27s = %d\n","Sidelenth of lattice",sidelength);
  else
  {
    printf("No Specified sidelength, defaulting to 10\n");
    sidelength = 10;
  }
  if(!config_lookup_bool(&cfg,"useclud",&useclud))
  {
    printf("No Specified update algorithm, defaulting to metropolis\n");
    useclud = 0;
  }
  if(useclud != 0)
    printf("%-27s = %s\n","Update Algorith","Cluster");
  else
    printf("%-27s = %s\n","Update Algorith","Metropolis");


  for(i = 0 ; i < 10 ; i++)
    printf("-");
  printf("\n");


  /* ALLOCATION OF LATTICE AND VECTORS */
  
  lattice = allocate_lattice(sidelength,spacedims,spindims);

  if(verbose_flag)
    fprintf(stderr,"Allocated %d points on the lattice\n",intpow(sidelength,spacedims));
  
  
  randomize_spins(lattice,sidelength,spacedims,spindims);

  if(verbose_flag) 
    print_lattice (lattice,sidelength,spacedims,spindims);





  /* CLEANUP */
  free_lattice(lattice,sidelength,spacedims);
  if(verbose_flag)
    fprintf(stderr,"Freed Lattice\n");


  return EXIT_SUCCESS;
}
