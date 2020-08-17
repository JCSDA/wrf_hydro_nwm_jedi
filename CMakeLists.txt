# (C) Copyright 2019-2019 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# WRF_HYDRO_NWM_JEDI
################################################################################

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

project( wrf_hydro_nwm_jedi C CXX Fortran)

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH})

set( CMAKE_DIRECOTYR_LABELS "wrf_hydro_nwm_jedi" )

set( ECBUILD_DEFAULT_BUILD_TYPE Release )
set( ENABLE_OS_TESTS           OFF CACHE BOOL "Disable OS tests" FORCE )
set( ENABLE_LARGE_FILE_SUPPORT OFF CACHE BOOL "Disable testing of large file support" FORCE )
#set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )

include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.5 )


################################################################################
# Project
################################################################################

ecbuild_declare_project()

ecbuild_enable_fortran( REQUIRED )
ecbuild_add_cxx11_flags()

set( WRF_HYDRO_NWM_JEDI_LINKER_LANGUAGE CXX )


################################################################################
# Dependencies
################################################################################
# Boost Headers
include_directories( ${Boost_INCLUDE_DIR} )

# NetCDF
find_package( NetCDF REQUIRED COMPONENTS Fortran )
include_directories( ${NETCDF_INCLUDE_DIRS} )

# eckit
ecbuild_use_package( PROJECT eckit VERSION 1.1.0 REQUIRED )
include_directories( ${ECKIT_INCLUDE_DIRS} )
link_libraries( ${ECKIT_LIBRARIES} )

# fckit
ecbuild_use_package( PROJECT fckit VERSION 0.4.1 REQUIRED )
include_directories( ${FCKIT_INCLUDE_DIRS} )

# oops
ecbuild_use_package( PROJECT oops VERSION 0.1.0 REQUIRED )
include_directories( ${OOPS_INCLUDE_DIRS} )

# ioda
ecbuild_use_package( PROJECT ioda VERSION 0.1.0 REQUIRED )
include_directories( ${IODA_INCLUDE_DIRS} )

# ufo
ecbuild_use_package( PROJECT ufo VERSION 0.1.0 REQUIRED )
include_directories( ${UFO_INCLUDE_DIRS} )


################################################################################
# Export package info
################################################################################
list( APPEND WRF_HYDRO_NWM_JEDI_TPLS LAPACK MPI NetCDF )

set( WRF_HYDRO_NWM_JEDI_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src
                       ${CMAKE_CURRENT_BINARY_DIR}/src
                       ${CMAKE_Fortran_MODULE_DIRECTORY} )
set( WRF_HYDRO_NWM_JEDI_LIBRARIES wrf_hydro_nwm_jedi )

get_directory_property( WRF_HYDRO_NWM_JEDI_DEFINITIONS COMPILE_DEFINITIONS )

foreach( _tpl ${WRF_HYDRO_NWM_JEDI_TPLS} )
string( TOUPPER ${_tpl} TPL )
list( APPEND WRF_HYDRO_NWM_JEDI_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
list( APPEND WRF_HYDRO_NWM_JEDI_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
list( APPEND WRF_HYDRO_NWM_JEDI_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()


################################################################################
# Sources
################################################################################
include( wrf_hydro_nwm-jedi_compiler_flags )
include_directories( ${WRF_HYDRO_NWM_JEDI_INCLUDE_DIRS} ${WRF_HYDRO_NWM_JEDI_EXTRA_INCLUDE_DIRS} )

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
