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

# ECMWF libs
# ----------
option("BUNDLE_SKIP_ECKIT" "Don't build eckit" "ON" ) # Skip eckit build unless user passes -DBUNDLE_SKIP_ECKIT=OFF
option("BUNDLE_SKIP_FCKIT" "Don't build fckit" "ON") # Build fckit unless user passes -DBUNDLE_SKIP_FCKIT=OFF
option("BUNDLE_SKIP_ATLAS" "Don't build atlas" "ON") # Build atlas unless user passes -DBUNDLE_SKIP_ATLAS=OFF

ecbuild_bundle( PROJECT eckit GIT "https://github.com/ecmwf/eckit.git" TAG 1.15.4 )
ecbuild_bundle( PROJECT fckit GIT "https://github.com/ecmwf/fckit.git" TAG 0.9.1 )
ecbuild_bundle( PROJECT atlas GIT "https://github.com/ecmwf/atlas.git" TAG 0.23.0 )

# External (required) observation operators
# ------------------------------
option("BUNDLE_SKIP_CRTM" "Don't build CRTM" "ON") # Skip crtm build unless user passes -DBUNDLE_SKIP_CRTM=OFF
ecbuild_bundle( PROJECT crtm GIT "https://github.com/jcsda-internal/crtm.git" BRANCH feature/ecbuild35 UPDATE )

# Core JEDI repositories
# ----------------------
ecbuild_bundle( PROJECT oops  GIT "https://github.com/jcsda-internal/oops.git"  BRANCH develop UPDATE )
ecbuild_bundle( PROJECT saber GIT "https://github.com/jcsda-internal/saber.git" BRANCH develop UPDATE )
ecbuild_bundle( PROJECT ioda  GIT "https://github.com/jcsda-internal/ioda.git"  BRANCH develop UPDATE )
ecbuild_bundle( PROJECT ufo   GIT "https://github.com/jcsda-internal/ufo.git"   BRANCH develop UPDATE )
ecbuild_bundle( PROJECT ufo-data  GIT "https://github.com/jcsda-internal/ufo-data.git"  BRANCH develop UPDATE )
ecbuild_bundle( PROJECT fv3-jedi-data  GIT "https://github.com/jcsda-internal/fv3-jedi-data.git"  BRANCH develop UPDATE )

# This repo
ecbuild_bundle( PROJECT wrf_hydro_nwm_jedi  SOURCE "../"  )

# Build Doxygen documentation
# ---------------------------
option(BUILD_FV3_BUNDLE_DOC "Build documentation" OFF)
if(BUILD_FV3_BUNDLE_DOC)
   add_subdirectory( Documentation )
endif(BUILD_FV3_BUNDLE_DOC)

ecbuild_bundle_finalize()

# cache repo hashes in a log file in the bundle.
execute_process ( COMMAND bash -c "${CMAKE_CURRENT_LIST_DIR}/log_repo_hashes.sh ${CMAKE_CURRENT_LIST_DIR}" )
