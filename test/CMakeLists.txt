list( APPEND wrf_hydro_nwm-jedi_test_input
  testinput/interface.yml
  )


# link the input files for the tests
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
foreach(FILENAME ${wrf_hydro_nwm-jedi_test_input})
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
    ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach()


#================================================================================
# Tests of class interfaces
#================================================================================

ecbuild_add_test(
  TARGET  test_wrf_hydro_nwm-jedi_geometry
  SOURCES executables/TestGeometry.cc
  ARGS    testinput/interface.yml
  LIBS    wrf_hydro_nwm-jedi )

# ecbuild_add_test(
#   TARGET  test_wrf_hydro_nwm-jedi_state
#   SOURCES executables/TestState.cc
#   ARGS    testinput/interface.yml
#   LIBS    wrf_hydro_nwm-jedi )

# ecbuild_add_test(
#   TARGET  test_wrf_hydro_nwm-jedi_modelauxcontrol
#   SOURCES executables/TestModelAuxControl.cc
#   ARGS    testinput/interface.yml
#   LIBS    wrf_hydro_nwm-jedi )

# ecbuild_add_test(
#   TARGET  test_wrf_hydro_nwm-jedi_model
#   SOURCES executables/TestModel.cc
#   ARGS    testinput/interface.yml
#   LIBS    wrf_hydro_nwm-jedi )

# ecbuild_add_test(
#   TARGET  test_wrf_hydro_nwm-jedi_increment
#   SOURCES executables/TestIncrement.cc
#   ARGS    testinput/interface.yml
#   LIBS    wrf_hydro_nwm-jedi )

# ecbuild_add_test(
#   TARGET  test_wrf_hydro_nwm-jedi_errorcovariance
#   SOURCES executables/TestErrorCovariance.cc
#   ARGS    testinput/interface.yml
#   LIBS    wrf_hydro_nwm-jedi )

# ecbuild_add_test(
#   TARGET  test_wrf_hydro_nwm-jedi_modelauxcovariance
#   SOURCES executables/TestModelAuxCovariance.cc
#   ARGS    testinput/interface.yml
#   LIBS    wrf_hydro_nwm-jedi )

# ecbuild_add_test(
#   TARGET  test_wrf_hydro_nwm-jedi_modelauxincrement
#   SOURCES executables/TestModelAuxIncrement.cc
#   ARGS    testinput/interface.yml
#   LIBS    wrf_hydro_nwm-jedi )

# ecbuild_add_test(
#   TARGET  test_wrf_hydro_nwm-jedi_linearmodel
#   SOURCES executables/TestLinearModel.cc
#   ARGS    testinput/interface.yml
#   LIBS    wrf_hydro_nwm-jedi )  

# ecbuild_add_test(
#   TARGET  test_wrf_hydro_nwm-jedi_geometryiterator
#   SOURCES executables/TestGeometryIterator.cc
#   ARGS    testinput/interface.yml
#   LIBS    wrf_hydro_nwm-jedi )  
