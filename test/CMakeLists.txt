list( APPEND wrf_hydro_nwm_jedi_test_input
  cp_file.sh
  testinput/geometry.yaml
  testinput/state.yaml
  testinput/interface.yml
  testinput/model.yaml
  testinput/increment.yaml
  testinput/error_covariance.yaml
  testinput/model_aux_control.yaml
  #  testinput/model_aux_covariance.yaml
  #  testinput/model_aux_increment.yaml
  testinput/linear_model.yaml
  testinput/make_obs_3d.yaml
  testinput/hofx3d.yaml
  testinput/staticb.yaml
  testinput/dirac.yaml
  testinput/3dvar.yaml
  testinput/makeobs.yaml
  testoutput/make_obs_3d.ref
  testoutput/hofx3d_identity.ref
  testoutput/staticb.ref
  testoutput/dirac.ref
  testoutput/3dvar.ref
  Data/wrf_hydro_nwm_files/wrfinput.nc
  Data/wrf_hydro_nwm_files/geometry_nwm_long_range_snow.nc
  Data/wrf_hydro_nwm_files/nwm_osse_ioda.nc
  Data/wrf_hydro_nwm_files/RESTART.2017010100_DOMAIN1
  Data/wrf_hydro_nwm_files/HYDRO_RST.2017-01-01_00:00_DOMAIN1
  Data/wrf_hydro_nwm_files/RESTART.2017010100_sneqv_plus5_DOMAIN1
  )

#list( APPEND wrf_hydro_nwm_jedi_test_input_cp
#  Data/wrf_hydro_nwm_files/RESTART.2017010100_sneqv_plus5_DOMAIN1
#  )


# link the input files for the tests
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/wrf_hydro_nwm_files)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/make_obs)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/hofx)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/staticb)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/dirac)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/3dvar)
foreach(FILENAME ${wrf_hydro_nwm_jedi_test_input})
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
    ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach()

#foreach(FILENAME ${wrf_hydro_nwm_jedi_test_input_cp})
#  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different
#    ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
#    ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
#endforeach()


#================================================================================
# Tests of class interfaces
#================================================================================

ecbuild_add_test(
  TARGET  test_wrf_hydro_nwm_jedi_geometry
  SOURCES executables/TestGeometry.cc
  ARGS    testinput/geometry.yaml
  LIBS    wrf_hydro_nwm_jedi )

ecbuild_add_test(
  TARGET  test_wrf_hydro_nwm_jedi_state
  SOURCES executables/TestState.cc
  ARGS    testinput/state.yaml
  LIBS    wrf_hydro_nwm_jedi )

ecbuild_add_test(
  TARGET  test_wrf_hydro_nwm_jedi_modelauxcontrol
  SOURCES executables/TestModelAuxControl.cc
  ARGS    testinput/model_aux_control.yaml
  LIBS    wrf_hydro_nwm_jedi )

ecbuild_add_test(
  TARGET  test_wrf_hydro_nwm_jedi_model
  SOURCES executables/TestModel.cc
  ARGS    testinput/model.yaml
  LIBS    wrf_hydro_nwm_jedi )

ecbuild_add_test(
  TARGET  test_wrf_hydro_nwm_jedi_increment
  SOURCES executables/TestIncrement.cc
  ARGS    testinput/increment.yaml
  LIBS    wrf_hydro_nwm_jedi )

ecbuild_add_test(
  TARGET  test_wrf_hydro_nwm_jedi_errorcovariance
  SOURCES executables/TestErrorCovariance.cc
  ARGS    testinput/error_covariance.yaml
  LIBS    wrf_hydro_nwm_jedi )

# ecbuild_add_test(
#   TARGET  test_wrf_hydro_nwm_jedi_modelauxcovariance
#   SOURCES executables/TestModelAuxCovariance.cc
#   ARGS    testinput/model_aux_covariance.yaml
#   LIBS    wrf_hydro_nwm_jedi )

# ecbuild_add_test(
#   TARGET  test_wrf_hydro_nwm_jedi_modelauxincrement
#   SOURCES executables/TestModelAuxIncrement.cc
#   ARGS    testinput/model_aux_increment.yaml
#   LIBS    wrf_hydro_nwm_jedi )

ecbuild_add_test(
  TARGET  test_wrf_hydro_nwm_jedi_linearmodel
  SOURCES executables/TestLinearModel.cc
  ARGS    testinput/linear_model.yaml
  LIBS    wrf_hydro_nwm_jedi )

# ecbuild_add_test(
#   TARGET  test_wrf_hydro_nwm_jedi_geometryiterator
#   SOURCES executables/TestGeometryIterator.cc
#   ARGS    testinput/interface.yml
#   LIBS    wrf_hydro_nwm_jedi )

ecbuild_add_test(
   TARGET       test_wrf_hydro_jedi_make_obs_3d
   TYPE         EXE
   COMMAND      ${CMAKE_BINARY_DIR}/bin/wrf_hydro_nwm_jedi_hofx.x
   DEPENDS      wrf_hydro_nwm_jedi_hofx.x
   ARGS         testinput/make_obs_3d.yaml
   MPI          1 )

ecbuild_add_test(
   TARGET       test_wrf_hydro_jedi_hofx3d_identity
   TYPE         EXE
   COMMAND      ${CMAKE_BINARY_DIR}/bin/wrf_hydro_nwm_jedi_hofx3d.x
   DEPENDS      wrf_hydro_nwm_jedi_hofx3d.x
   ARGS         testinput/hofx3d.yaml
   MPI          1 )

ecbuild_add_test(
   TARGET       test_wrf_hydro_jedi_stage_lsm
   TYPE         SCRIPT
   COMMAND      ${BASH_PROGRAM}
   ARGS         ${CMAKE_CURRENT_BINARY_DIR}/cp_file.sh
                Data/wrf_hydro_nwm_files/RESTART.2017010100_sneqv_plus5_DOMAIN1
                Data/staticb/stddev.lsm.an.2017-01-01T00:00:00Z
                Data/staticb/rh.lsm.an.2017-01-01T00:00:00Z
                Data/dirac/dirac_SABER.lsm.an.2017-01-01T00:00:00Z
                Data/3dvar/3dvar.lsm.an.2017-01-01T00:00:00Z )

ecbuild_add_test(
   TARGET       test_wrf_hydro_jedi_stage_hydro
   TYPE         SCRIPT
   COMMAND      ${BASH_PROGRAM}
   ARGS         ${CMAKE_CURRENT_BINARY_DIR}/cp_file.sh
                Data/wrf_hydro_nwm_files/HYDRO_RST.2017-01-01_00:00_DOMAIN1
                Data/staticb/stddev.hyd.an.2017-01-01T00:00:00Z
                Data/staticb/rh.hyd.an.2017-01-01T00:00:00Z
                Data/dirac/dirac_SABER.hyd.an.2017-01-01T00:00:00Z
                Data/3dvar/3dvar.hyd.an.2017-01-01T00:00:00Z )

ecbuild_add_test(
   TARGET       test_wrf_hydro_jedi_staticb
   TYPE         EXE
   COMMAND      ${CMAKE_BINARY_DIR}/bin/wrf_hydro_nwm_jedi_error_covariance_training.x
   DEPENDS      wrf_hydro_nwm_jedi_error_covariance_training.x
   ARGS         testinput/staticb.yaml
   MPI          1 )

ecbuild_add_test(
   TARGET       test_wrf_hydro_jedi_dirac
   TYPE         EXE
   COMMAND      ${CMAKE_BINARY_DIR}/bin/wrf_hydro_nwm_jedi_dirac.x
   DEPENDS      wrf_hydro_nwm_jedi_dirac.x
   ARGS         testinput/dirac.yaml
   MPI          1 )

ecbuild_add_test(
   TARGET       test_wrf_hydro_jedi_var
   TYPE         EXE
   COMMAND      ${CMAKE_BINARY_DIR}/bin/wrf_hydro_nwm_jedi_var.x
   DEPENDS      wrf_hydro_nwm_jedi_var.x
   ARGS         testinput/3dvar.yaml
   MPI          1 )
