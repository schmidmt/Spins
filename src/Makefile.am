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


bin_PROGRAMS = spinsMetropolisScanT spinsClusterScanT spinsClusterSingleBeta spinsCriticalExps

spinsMetropolisScanT_SOURCES = spinsMetropolisScanT.c clusterupdate.c metropolis.c lattice.c common.c physics.c
spinsMetropolisScanT_CFLAGS = -pg -Wall  -DHAVE_INLINE

spinsClusterScanT_SOURCES = spinsClusterScanT.c clusterupdate.c lattice.c common.c physics.c
spinsClusterScanT_CFLAGS = -pg -Wall  -DHAVE_INLINE

spinsClusterSingleBeta_SOURCES = spinsClusterSingleBeta.c clusterupdate.c lattice.c common.c physics.c
spinsClusterSingleBeta_CFLAGS = -pg -Wall  -DHAVE_INLINE

spinsCriticalExps_SOURCES = spinsCriticalExps.c clusterupdate.c lattice.c common.c physics.c
spinsCriticalExps_CFLAGS = -pg -Wall  -DHAVE_INLINe


AM_CPPFLAGS = -I$(top_srcdir)/includes

EXTRA_DIST = example2.cfg
