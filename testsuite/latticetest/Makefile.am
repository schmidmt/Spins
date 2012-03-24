#    Copyright 2012 Michael Schmidt (mts@colorado.edu)
#    #
#    This file is part of spins.
#    #
#    spins is free software: you can redistribute it and/or modify
#    #    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    #    (at your option) any later version.
#
##    spins is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    #
#    You should have received a copy of the GNU General Public License
#    #    along with spins.  If not, see <http://www.gnu.org/licenses/>.

AUTOMAKE_OPTIONS = dejagnu
bin_PROGRAMS = latticetest
latticetest_SOURCES = latticetest.c $(top_srcdir)/src/lattice.c $(top_srcdir)/src/common.c

RUNTESTFLAGS = --tool latticetest LATTICETEST=`pwd`/latticetest --srcdir $(top_srcdir)/testsuite

AM_CPPFLAGS = -I$(top_srcdir)/includes
