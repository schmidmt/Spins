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


/**************
 * Parse argv *
 **************/
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
	strcpy (conf_file, optarg);
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


if (!specified_file_flag)
  {
    fprintf (stderr,
	     "#No file specified, please specify one with the -f flag\n");
    exit (EXIT_FAILURE);
  }


/***********************
 * Reading Config File *
 ************************/
if (!config_read_file (&cfg, conf_file))
  {
    fprintf (stderr, "%s:%d - %s\n", config_error_file (&cfg),
	     config_error_line (&cfg), config_error_text (&cfg));
    config_destroy (&cfg);
    return (EXIT_FAILURE);
  }

if (config_lookup_string (&cfg, "name", &str))
  {
    strcpy (outputfile_name, str);
    strcat (outputfile_name, ".dat");
    if ((outputfp = fopen (outputfile_name, "w+")) == NULL)
      {
	fprintf (stderr, "Cannot open file %s\n", outputfile_name);
	exit (EXIT_FAILURE);
      }
    for (i = 0; i < 10; i++)
      {
	printf ("#");
	fprintf (outputfp, "#");
      }
    printf ("\n");
    fprintf (outputfp, "\n");
    printf ("#%-27s = '%s'\n", "Run Name", str);
    fprintf (outputfp, "#%-27s = '%s'\n", "Run Name", str);
  }
else
  {
    fprintf (stderr, "#No 'name' setting found in configuration file.\n");
    return (EXIT_FAILURE);
  }

if (config_lookup_int (&cfg, "spindims", &conf.spindims))
  {
    printf ("#%-27s = %d\n", "Number of spin dimensions", conf.spindims);
    fprintf (outputfp, "#%-27s = %d\n", "Number of spin dimensions",
	     conf.spindims);
  }
else
  {
    printf ("#No Specified spindims, defaulting to 2\n");
    conf.spindims = 2;
  }

if (config_lookup_int (&cfg, "spacedims", &conf.spacedims))
  {
    printf ("#%-27s = %d\n", "Number of space dimensions", conf.spacedims);
    fprintf (outputfp, "#%-27s = %d\n", "Number of space dimensions",
	     conf.spacedims);
  }
else
  {
    printf ("#No Specified spacedims, defaulting to 2\n");
    conf.spacedims = 2;
  }

if (config_lookup_int (&cfg, "sidelength", &conf.sidelength))
  {
    printf ("#%-27s = %d\n", "Sidelenth of lattice", conf.sidelength);
    fprintf (outputfp, "#%-27s = %d\n", "Sidelenth of lattice",
	     conf.sidelength);
  }
else
  {
    printf ("#No Specified sidelength, defaulting to 10\n");
    conf.sidelength = 10;
  }

if (!config_lookup_bool (&cfg, "randomize_spins", &random_spins))
  {
    printf ("#Randomize_spins not specified, defaulting to random.\n");
    random_spins = 0;
  }

if (random_spins == 1)
  {
    printf ("#%-27s = %s\n", "Random Spins", "True");
    fprintf (outputfp, "#%-27s = %s\n", "Random Spins", "True");
  }
else
  {
    printf ("#%-27s = %s\n", "Random Spins", "False");
    fprintf (outputfp, "#%-27s = %s\n", "Random Spins", "False");
  }

if (!config_lookup_bool (&cfg, "useclud", &useclud))
  {
    printf ("#No Specified update algorithm, defaulting to metropolis\n");
    useclud = 0;
  }

if (useclud != 0)
  {
    printf ("#%-27s = %s\n", "Update Algorith", "Cluster");
    fprintf (outputfp, "#%-27s = %s\n", "Update Algorith", "Cluster");
  }
else
  {
    printf ("#%-27s = %s\n", "Update Algorith", "Metropolis");
    fprintf (outputfp, "#%-27s = %s\n", "Update Algorith", "Metropolis");
  }

if (config_lookup_float (&cfg, "beta_start", &beta_start))
  {
    printf ("#%-27s = %e\n", "Beta Start", beta_start);
    fprintf (outputfp, "#%-27s = %e\n", "Beta Start", beta_start);
  }
else
  {
    printf ("#No Specified beta_start, defaulting to 0\n");
    beta_start = 1;
  }

if (config_lookup_float (&cfg, "beta_end", &beta_end))
  {
    printf ("#%-27s = %e\n", "Beta End", beta_end);
    fprintf (outputfp, "#%-27s = %e\n", "Beta End", beta_end);
  }
else
  {
    printf ("#No Specified beta_end, defaulting to 10\n");
    beta_end = 10;
  }

if (config_lookup_float (&cfg, "beta_step", &beta_step))
  {
    printf ("#%-27s = %e\n", "Beta Step", beta_step);
    fprintf (outputfp, "#%-27s = %e\n", "Beta Step", beta_step);
  }
else
  {
    printf ("#No Specified beta_step, defaulting to 0.1\n");
    beta_step = 0.1;
  }

if (config_lookup_int (&cfg, "max_settle", &conf.max_settle))
  {
    printf ("#%-27s = %d\n", "Max Settle Iterations", conf.max_settle);
    fprintf (outputfp, "#%-27s = %d\n", "Max Settle Iterations",
	     conf.max_settle);
  }
else
  {
    printf ("#No Specified max_settle, defaulting to 1000\n");
    conf.max_settle = 10000;
  }

if (config_lookup_int (&cfg, "rng_seed", &rng_seed))
  {
    printf ("#%-27s = %d\n", "Random Seed", rng_seed);
    fprintf (outputfp, "#%-27s = %d\n", "Random Seed", rng_seed);
    gsl_rng_set (conf.rng, rng_seed);
  }
if (config_lookup_int (&cfg, "block_size", &conf.block_size))
  {
    printf ("#%-27s = %d\n", "Blocking Size", conf.block_size);
    fprintf (outputfp, "#%-27s = %d\n", "Blocking Size", conf.block_size);
  }
else
  {
    printf ("#No Specified block_size, defaulting to 1000\n");
    conf.block_size = 1000;
  }

if (config_lookup_int (&cfg, "blocks", &conf.blocks))
  {
    printf ("#%-27s = %d\n", "Blocks", conf.blocks);
    fprintf (outputfp, "#%-27s = %d\n", "Blocks", conf.blocks);
  }
else
  {
    printf ("#No Specified blocks, defaulting to 10\n");
    conf.blocks = 10;
  }

for (i = 0; i < 10; i++)
  {
    printf ("#");
    fprintf (outputfp, "#");
  }

printf ("\n");
fprintf (outputfp, "\n");
