AUTOMAKE_OPTIONS = dejagnu
bin_PROGRAMS = physicstest
physicstest_SOURCES = physicstest.c $(top_srcdir)/src/physics.c $(top_srcdir)/src/common.c $(top_srcdir)/src/lattice.c

RUNTESTFLAGS = --tool physicstest TEST=`pwd`/physicstest --srcdir $(top_srcdir)/testsuite

AM_CPPFLAGS = -I$(top_srcdir)/includes
