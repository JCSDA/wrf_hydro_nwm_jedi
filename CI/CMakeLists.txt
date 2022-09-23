# (C) COPYRIGHT 2018-2019 UCAR:
#
# THIS SOFTWARE IS LICENSED UNDER THE TERMS OF THE APACHE LICENCE VERSION 2.0
# WHICH CAN BE OBTAINED AT HTTP://WWW.APACHE.ORG/LICENSES/LICENSE-2.0.

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

find_package( ecbuild 3.5 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild)

# WRF_HYDRO_NWM_JEDI
project( WRF_HYDRO_NWM_JEDI VERSION 1.0.0 LANGUAGES C CXX Fortran )

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include( ecbuild_bundle )

# set( ECBUILD_DEFAULT_BUILD_TYPE Release )
set( ECBUILD_DEFAULT_BUILD_TYPE RelWithDebInfo )

# Enable MPI
set( ENABLE_MPI ON CACHE BOOL "Compile with MPI")

# Initialize
ecbuild_bundle_initialize()

# Add git functions from jedi-cmake
include( $ENV{jedi_cmake_ROOT}/share/jedicmake/Functions/git_functions.cmake )

# External (required) observation operators
# ------------------------------
option("BUNDLE_SKIP_CRTM" "Don't build CRTM" "ON") # Skip crtm build unless user passes -DBUNDLE_SKIP_CRTM=OFF
ecbuild_bundle( PROJECT crtm GIT "https://github.com/jcsda-internal/crtm.git" )

# Core JEDI repositories
# ----------------------
ecbuild_bundle( PROJECT oops  GIT "https://github.com/jcsda-internal/oops.git" )
ecbuild_bundle( PROJECT saber GIT "https://github.com/jcsda-internal/saber.git" )
ecbuild_bundle( PROJECT ioda  GIT "https://github.com/jcsda-internal/ioda.git" )
ecbuild_bundle( PROJECT ufo   GIT "https://github.com/jcsda-internal/ufo.git" )
ecbuild_bundle( PROJECT ufo-data  GIT "https://github.com/jcsda-internal/ufo-data.git" )

# This repo
ecbuild_bundle( PROJECT wrf_hydro_nwm_jedi  GIT "https://github.com/jcsda-internal/wrf_hydro_nwm_jedi.git" )

# Build Doxygen documentation
# ---------------------------
option(BUILD_FV3_BUNDLE_DOC "Build documentation" OFF)
if(BUILD_FV3_BUNDLE_DOC)
   add_subdirectory( Documentation )
endif(BUILD_FV3_BUNDLE_DOC)

ecbuild_bundle_finalize()

include(cmake/cdash-integration.cmake)
include(CTest)

# cache repo hashes in a log file in the bundle.
find_program( BASH_PROGRAM bash )
execute_process ( COMMAND ${BASH_PROGRAM} -c "${CMAKE_CURRENT_LIST_DIR}/log_repo_hashes.sh ${CMAKE_CURRENT_LIST_DIR}" )
