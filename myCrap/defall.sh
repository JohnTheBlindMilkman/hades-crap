#!/bin/sh
###############################################################################
#
#  Hades Software Environment Setup Script (Bourne shell family version)
#
#  Each Software Collection should have it dedicated verion of this script.
#  They can be distinguished e.g. by the used paths.
#
#  Author: Simon Lang, GSI, 30.09.2006
#
###############################################################################

. /cvmfs/hadessoft.gsi.de/install/debian10/admin/hsc-functions.sh

hsc_checkContext

# Root, and dependent packages
export ROOTSYS=/cvmfs/hadessoft.gsi.de/install/debian10/root-6.24.02

# Global Hydra Location
export HADDIR=/cvmfs/hadessoft.gsi.de/install/debian10/6.24.02/hydra2-6.2
export ROOT_INCLUDE_PATH=${HADDIR}/include

# Private Hydra Location - not used by default
#export MYHADDIR=/e.g./somewhere/in/my/home/directory
#export ROOT_INCLUDE_PATH=${MYHADDIR}/include:${HADDIR}/include
export MYHADDIR=/misc/hadesprojects/analysis_workshop/2022/hydra/own_lib
export ROOT_INCLUDE_PATH=${MYHADDIR}/include:${HADDIR}/include

# Oracle
export ORACLE_HOME=/cvmfs/hadessoft.gsi.de/install/debian10/oracle/client
export ORA_USER=hades_ana/hades@db-hades

# CERNLIB - for HGeant
export CERN_ROOT=/cvmfs/hadessoft.gsi.de/install/debian10/cernlib_gfortran/2005

# PLUTO
export PLUTODIR=

# RFIO support
export ADSM_BASE_NEW=/misc/hadessoftware/etch32/install/gstore-may07

hsc_setEnvironment
hsc_shortPrintout

hsc_finalizeScript
