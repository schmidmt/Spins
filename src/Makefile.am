bin_PROGRAMS = spinsMetropolisScanT spinsClusterScanT spinsClusterSingleBeta

spinsMetropolisScanT_SOURCES = spinsMetropolisScanT.c clusterupdate.c metropolis.c lattice.c common.c physics.c
spinsMetropolisScanT_CFLAGS = -pg -Wall  -DHAVE_INLINE

spinsClusterScanT_SOURCES = spinsClusterScanT.c clusterupdate.c lattice.c common.c physics.c
spinsClusterScanT_CFLAGS = -pg -Wall  -DHAVE_INLINE

spinsClusterSingleBeta_SOURCES = spinsClusterSingleBeta.c clusterupdate.c lattice.c common.c physics.c
spinsClusterSingleBeta_CFLAGS = -pg -Wall  -DHAVE_INLINE


AM_CPPFLAGS = -I$(top_srcdir)/includes

EXTRA_DIST = example2.cfg
