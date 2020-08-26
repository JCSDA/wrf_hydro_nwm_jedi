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


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Optional repositories
if ( BUILD_ECKIT )
   ecbuild_bundle( PROJECT eckit         GIT "https://github.com/JCSDA/eckit.git"           UPDATE BRANCH release-stable )
endif ()


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# required repositories: sets for different branches when things get out of sync

# -----------------------------------------------------------------------------
# Bleeding-edge
# ecbuild_bundle( PROJECT fckit           GIT "https://github.com/JCSDA/fckit.git"           UPDATE BRANCH release-stable )
# ecbuild_bundle( PROJECT atlas           GIT "https://github.com/JCSDA/atlas.git"           UPDATE BRANCH release-stable )
# ecbuild_bundle( PROJECT oops            GIT "https://github.com/JCSDA/oops.git"            UPDATE BRANCH develop )
# ecbuild_bundle( PROJECT saber           GIT "https://github.com/JCSDA/saber.git"           UPDATE BRANCH develop )
# ecbuild_bundle( PROJECT ioda            GIT "https://github.com/JCSDA/ioda.git"            UPDATE BRANCH develop )
# ecbuild_bundle( PROJECT ioda-converters GIT "https://github.com/JCSDA/ioda-converters.git" UPDATE BRANCH develop )
# ecbuild_bundle( PROJECT ufo             GIT "https://github.com/JCSDA/ufo.git"             UPDATE BRANCH develop )

# -----------------------------------------------------------------------------
# obs_trait: July/August 2020
ecbuild_bundle( PROJECT fckit           GIT "https://github.com/JCSDA/fckit.git"           UPDATE BRANCH release-stable )
ecbuild_bundle( PROJECT atlas           GIT "https://github.com/JCSDA/atlas.git"           UPDATE BRANCH release-stable )
ecbuild_bundle( PROJECT oops            GIT "https://github.com/JCSDA/oops.git"            UPDATE TAG    24f7b4815f6dda0eb9116611570f69025f202d72 )
ecbuild_bundle( PROJECT saber           GIT "https://github.com/JCSDA/saber.git"           UPDATE TAG    e0998b6bdfb6f63b111ee6cf60e31b37189fea3a )
ecbuild_bundle( PROJECT ioda            GIT "https://github.com/JCSDA/ioda.git"            UPDATE TAG    738fa41f0335226941f50fa57bed941be5cf520c )
ecbuild_bundle( PROJECT ioda-converters GIT "https://github.com/JCSDA/ioda-converters.git" UPDATE BRANCH develop )
ecbuild_bundle( PROJECT ufo             GIT "https://github.com/JCSDA/ufo.git"             UPDATE TAG    f894165310056daa2ca584570b9274a6195b11cd )

# -----------------------------------------------------------------------------
# feature/3dvar: July/August 2020.
# singularity container: jedi-gnu-openmpi-dev_sha256.de5e28890fee511bf22db77b0083c1255b5b96c572ff1ce05e7f34c4c99d76cf.sif
# ecbuild_bundle( PROJECT fckit           GIT "https://github.com/JCSDA/fckit.git"           UPDATE TAG 2fad06b133744751b91cf9d0fb5d5b7de7ff1c57 )
# ecbuild_bundle( PROJECT atlas           GIT "https://github.com/JCSDA/atlas.git"           UPDATE TAG d075f57b7cf561e082a2904261801447758bd5cf )
# ecbuild_bundle( PROJECT oops            GIT "https://github.com/JCSDA/oops.git"            UPDATE TAG 0566cc7ed432f442e3d912953f77191f0b605bb2 )
# ecbuild_bundle( PROJECT saber           GIT "https://github.com/JCSDA/saber.git"           UPDATE TAG e9518dcce4115c8d07658e8c4400cbbc87bc5e5a )
# ecbuild_bundle( PROJECT ioda            GIT "https://github.com/JCSDA/ioda.git"            UPDATE TAG 3e49e5322b04eebba32de3b1efdb7fc971d2dc43 )
# ecbuild_bundle( PROJECT ioda-converters GIT "https://github.com/JCSDA/ioda-converters.git" UPDATE TAG f57c8b0d6b2ff974438b1aecf3fa97a5d667f602 )
# ecbuild_bundle( PROJECT ufo             GIT "https://github.com/JCSDA/ufo.git"             UPDATE TAG cddd7c949a64e0430ebb3e1ac9010f95155fda77 )


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# This repo
ecbuild_bundle( PROJECT wrf_hydro_nwm_jedi      SOURCE "../"             )


ecbuild_bundle_finalize()


# cache repo hashes in a log file in the bundle.
execute_process ( COMMAND bash -c "${CMAKE_CURRENT_LIST_DIR}/log_repo_hashes.sh ${CMAKE_CURRENT_LIST_DIR}" )
