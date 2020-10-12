list( APPEND wrf_hydro_nwm_jedi_test_input
  compare.sh
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
  testinput/hofx_nomodel.yaml
  testinput/3dvar.yaml
  testinput/makeobs.yaml
  testoutput/hofx_nomodel_identity.ref
  testoutput/3dvar.ref
  Data/wrf_hydro_nwm_files/wrfinput.nc
  Data/wrf_hydro_nwm_files/geometry_nwm_long_range_snow.nc
  Data/wrf_hydro_nwm_files/RESTART.2017010100_DOMAIN1
  Data/wrf_hydro_nwm_files/HYDRO_RST.2017-01-01_00:00_DOMAIN1
  Data/wrf_hydro_nwm_files/nwm_osse_ioda.nc
  )


# link the input files for the tests
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/wrf_hydro_nwm_files)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/make_obs)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/hofx)
foreach(FILENAME ${wrf_hydro_nwm_jedi_test_input})
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
    ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach()


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
   ARGS         testinput/make_obs_3d.yaml testoutput/make_obs_3d.log.out
   MPI          1 )

ecbuild_add_test(
   TARGET       test_wrf_hydro_jedi_hofx_nomodel_identity_exe
   TYPE         EXE
   COMMAND      ${CMAKE_BINARY_DIR}/bin/wrf_hydro_nwm_jedi_hofxnomodel.x
   ARGS         testinput/hofx_nomodel.yaml testoutput/hofx_nomodel_identity.log.out
   MPI          1 )

ecbuild_add_test(
   TARGET       test_wrf_hydro_jedi_hofx_nomodel_identity_cmp
   TYPE		SCRIPT
   COMMAND	"${CMAKE_CURRENT_BINARY_DIR}/compare.sh"
   ARGS         testoutput/hofx_nomodel_identity.log.out testoutput/hofx_nomodel_identity.ref)

ecbuild_add_test(
   TARGET       test_wrf_hydro_jedi_var_exe
   TYPE         EXE
   COMMAND      ${CMAKE_BINARY_DIR}/bin/wrf_hydro_nwm_jedi_var.x
   ARGS         testinput/3dvar.yaml testoutput/3dvar.log.out
   MPI          1 )

ecbuild_add_test(
   TARGET       test_wrf_hydro_jedi_var_cmp
   TYPE		SCRIPT
   COMMAND	"${CMAKE_CURRENT_BINARY_DIR}/compare.sh"
   ARGS         testoutput/3dvar.log.out testoutput/3dvar.ref)
