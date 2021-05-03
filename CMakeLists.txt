# (C) Copyright 2019-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# WRF_HYDRO_NWM_JEDI
################################################################################

cmake_minimum_required( VERSION 3.12 )
project( wrf_hydro_nwm_jedi VERSION 0.1 LANGUAGES C CXX Fortran)

# ecbuild integratoin
find_package(ecbuild 3.3.2 REQUIRED)
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()

set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include( wrf_hydro_nwm_jedi_compiler_flags )

################################################################################
# Dependencies
################################################################################
# Boost Headers
include_directories( ${Boost_INCLUDE_DIR} )

find_package( NetCDF REQUIRED COMPONENTS Fortran )
find_package( eckit  1.11.6     REQUIRED )
find_package( fckit  0.7.0      REQUIRED )
find_package( atlas  0.20.2     REQUIRED )
find_package( oops   0.1.0      REQUIRED )
find_package( saber  0.0.1      REQUIRED )
find_package( ioda   0.1.0      REQUIRED )
find_package( ufo    0.1.0      REQUIRED )

################################################################################
# Sources
################################################################################

set( WRF_HYDRO_NWM_JEDI_LINKER_LANGUAGE CXX )

add_subdirectory( src )
add_subdirectory( test )

# Build Doxygen documentation
# add_subdirectory( doxygen )

if(ECBUILD_INSTALL_FORTRAN_MODULES)
  install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/${CMAKE_CFG_INTDIR} DESTINATION ${INSTALL_INCLUDE_DIR} )
endif()

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME wrf_hydro_nwm_jedi )

# print the summary of the configuration
ecbuild_print_summary()
