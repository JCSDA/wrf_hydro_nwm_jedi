# (C) COPYRIGHT 2018-2019 UCAR
#
# THIS SOFTWARE IS LICENSED UNDER THE TERMS OF THE APACHE LICENCE VERSION 2.0
# WHICH CAN BE OBTAINED AT HTTP://WWW.APACHE.ORG/LICENSES/LICENSE-2.0.

# WRF_HYDRO_NWM_JEDI
project( WRF_HYDRO_NWM_JEDI C CXX Fortran )

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )
include( ecbuild_bundle )

set( ECBUILD_DEFAULT_BUILD_TYPE RelWithDebInfo )
# set( ECBUILD_DEFAULT_BUILD_TYPE Release )
set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

ecbuild_bundle_initialize()
ecbuild_requires_macro_version( 2.7 )

#ecbuild_add_option( FEATURE OMP
#                    DEFAULT ON
#                    DESCRIPTION "Use OpenMP" )

ecbuild_add_option( FEATURE BUILD_ECKIT
                    DEFAULT OFF
                    DESCRIPTION "download and build eckit (not needed if in a jedi container)" )


# Optional repositories

if ( BUILD_ECKIT )
   ecbuild_bundle( PROJECT eckit         GIT "https://github.com/JCSDA/eckit.git"           UPDATE BRANCH release-stable )
endif ()


# required repositories

ecbuild_bundle( PROJECT fckit           GIT "https://github.com/JCSDA/fckit.git"           UPDATE BRANCH release-stable )
ecbuild_bundle( PROJECT atlas           GIT "https://github.com/JCSDA/atlas.git"           UPDATE BRANCH release-stable )
ecbuild_bundle( PROJECT saber           GIT "https://github.com/JCSDA/saber.git"           UPDATE BRANCH develop )
ecbuild_bundle( PROJECT oops            GIT "https://github.com/JCSDA/oops.git"            UPDATE BRANCH develop )
ecbuild_bundle( PROJECT ioda            GIT "https://github.com/JCSDA/ioda.git"            UPDATE BRANCH develop )
ecbuild_bundle( PROJECT ioda-converters GIT "https://github.com/JCSDA/ioda-converters.git" UPDATE BRANCH develop )
ecbuild_bundle( PROJECT ufo             GIT "https://github.com/JCSDA/ufo.git"             UPDATE BRANCH develop )
ecbuild_bundle( PROJECT wrf_hydro_nwm_jedi      SOURCE "../"             )
  
ecbuild_bundle_finalize()
