AC_INIT([spins],[0.1],[mts@colorado.edu])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

AC_CHECK_LIB([config],[config_read_file])

AC_PROG_CC
AM_PROG_CC_C_O
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile testsuite/Makefile testsuite/latticetest/Makefile testsuite/physicstest/Makefile])
AC_OUTPUT
