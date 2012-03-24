#    Copyright 2012 Michael Schmidt (mts@colorado.edu)
#
#    This file is part of spins.
#
#    spins is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    spins is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with spins.  If not, see <http://www.gnu.org/licenses/>.

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
