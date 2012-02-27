AUTOMAKE_OPTIONS = dejagnu
bin_PROGRAMS = latticetest
latticetest_SOURCES = latticetest.c $(top_srcdir)/src/lattice.c $(top_srcdir)/src/common.c

RUNTESTFLAGS = --tool latticetest LATTICETEST=`pwd`/latticetest --srcdir $(top_srcdir)/testsuite

AM_CPPFLAGS = -I$(top_srcdir)/includes
