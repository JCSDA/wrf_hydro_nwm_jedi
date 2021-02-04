# (C) COPYRIGHT 2018-2019 UCAR:
#
# THIS SOFTWARE IS LICENSED UNDER THE TERMS OF THE APACHE LICENCE VERSION 2.0
# WHICH CAN BE OBTAINED AT HTTP://WWW.APACHE.ORG/LICENSES/LICENSE-2.0.

# WRF_HYDRO_NWM_JEDI
project( WRF_HYDRO_NWM_JEDI C CXX Fortran )

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

find_package(ecbuild)

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


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Optional repositories
if ( BUILD_ECKIT )
   ecbuild_bundle( PROJECT eckit         GIT "https://github.com/JCSDA/eckit.git"           UPDATE BRANCH release-stable )
endif ()

# required repositories: sets for different branches when things get out of sync
# Formula for specific commits. See recent hashes in .repo_hash_log
# ecbuild_bundle( PROJECT repo            GIT "https://github.com/JCSDA-internal/repo.git"          UPDATE TAG    24f7b4815f6dda0eb9116611570f69025f202d72 )
ecbuild_bundle( PROJECT fckit           GIT "https://github.com/JCSDA-internal/fckit.git" )
ecbuild_bundle( PROJECT atlas           GIT "https://github.com/JCSDA-internal/atlas.git" )
ecbuild_bundle( PROJECT ioda-converters GIT "https://github.com/JCSDA-internal/ioda-converters.git" )
ecbuild_bundle( PROJECT oops            GIT "https://github.com/JCSDA-internal/oops.git" )
ecbuild_bundle( PROJECT saber           GIT "https://github.com/JCSDA-internal/saber.git" )
ecbuild_bundle( PROJECT ioda            GIT "https://github.com/JCSDA-internal/ioda.git" )
ecbuild_bundle( PROJECT ufo             GIT "https://github.com/JCSDA-internal/ufo.git" )

# This repo
ecbuild_bundle( PROJECT wrf_hydro_nwm_jedi  GIT "https://github.com/JCSDA-internal/wrf_hydro_nwm_jedi" )


ecbuild_bundle_finalize()


# cache repo hashes in a log file in the bundle.
execute_process ( COMMAND bash -c "${CMAKE_CURRENT_LIST_DIR}/log_repo_hashes.sh ${CMAKE_CURRENT_LIST_DIR}" )
