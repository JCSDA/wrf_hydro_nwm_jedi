list( APPEND wrf_hydro_nwm_jedi_test_input
  testinput/interface.yml
  testinput/hofx_nomodel.yaml
  testoutput/hofx_nomodel_identity.ref
  Data/wrf_hydro_nwm_files/wrfinput.nc
  Data/wrf_hydro_nwm_files/RESTART.2017010100_DOMAIN1
  Data/wrf_hydro_nwm_files/nwm_osse_ioda.nc
  )


# link the input files for the tests
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testoutput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Data/wrf_hydro_nwm_files)
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
  ARGS    testinput/interface.yml
  LIBS    wrf_hydro_nwm_jedi )

ecbuild_add_test(
  TARGET  test_wrf_hydro_nwm_jedi_state
  SOURCES executables/TestState.cc
  ARGS    testinput/interface.yml
  LIBS    wrf_hydro_nwm_jedi )

ecbuild_add_test(
  TARGET  test_wrf_hydro_nwm_jedi_modelauxcontrol
  SOURCES executables/TestModelAuxControl.cc
  ARGS    testinput/interface.yml
  LIBS    wrf_hydro_nwm_jedi )

# ecbuild_add_test(
#   TARGET  test_wrf_hydro_nwm_jedi_model
#   SOURCES executables/TestModel.cc
#   ARGS    testinput/interface.yml
#   LIBS    wrf_hydro_nwm_jedi )

ecbuild_add_test(
  TARGET  test_wrf_hydro_nwm_jedi_increment
  SOURCES executables/TestIncrement.cc
  ARGS    testinput/interface.yml
  LIBS    wrf_hydro_nwm_jedi )

ecbuild_add_test(
  TARGET  test_wrf_hydro_nwm_jedi_errorcovariance
  SOURCES executables/TestErrorCovariance.cc
  ARGS    testinput/interface.yml
  LIBS    wrf_hydro_nwm_jedi )

ecbuild_add_test(
  TARGET  test_wrf_hydro_nwm_jedi_modelauxcovariance
  SOURCES executables/TestModelAuxCovariance.cc
  ARGS    testinput/interface.yml
  LIBS    wrf_hydro_nwm_jedi )

ecbuild_add_test(
  TARGET  test_wrf_hydro_nwm_jedi_modelauxincrement
  SOURCES executables/TestModelAuxIncrement.cc
  ARGS    testinput/interface.yml
  LIBS    wrf_hydro_nwm_jedi )

ecbuild_add_test(
  TARGET  test_wrf_hydro_nwm_jedi_linearmodel
  SOURCES executables/TestLinearModel.cc
  ARGS    testinput/interface.yml
  LIBS    wrf_hydro_nwm_jedi )

# ecbuild_add_test(
#   TARGET  test_wrf_hydro_nwm_jedi_geometryiterator
#   SOURCES executables/TestGeometryIterator.cc
#   ARGS    testinput/interface.yml
#   LIBS    wrf_hydro_nwm_jedi )

ecbuild_add_test(
   TARGET       test_wrf_hydro_jedi_hofx_nomodel_identity_exe
   TYPE         EXE
   COMMAND      ${CMAKE_BINARY_DIR}/bin/wrf_hydro_nwm_jedi_hofxnomodel.x
   ARGS         testinput/hofx_nomodel.yaml testoutput/hofx_nomodel_identity.run
   MPI          1 )

ecbuild_add_test(
   TARGET       test_wrf_hydro_jedi_hofx_nomodel_identity_cmp
   TYPE		SCRIPT
   COMMAND	"${CMAKE_BINARY_DIR}/bin/compare.sh"
   ARGS         testoutput/hofx_nomodel_identity.run testoutput/hofx_nomodel_identity.ref)
#   COMMAND      ${CMAKE_BINARY_DIR}/bin/wrf_hydro_nwm_jedi_hofxnomodel.x
#   ARGS         testinput/hofx_nomodel.yaml testoutput/hofx_nomodel_identity.run
#   MPI          1 )
