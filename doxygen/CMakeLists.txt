# first we can indicate the documentation build as an option and
# set it to OFF by default since documenation will usually be
# built at the bundle level
option(BUILD_WRF_HYDRO_NWM_JEDI_DOC "Build WRF_HYDRO_NWM_JEDI documentation" OFF)

if (BUILD_WRF_HYDRO_NWM_JEDI_DOC)
    message("WRF_HYDRO_NWM_JEDI Documentation Requested")
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
        add_custom_target( doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating WRF_HYDRO_NWM_JEDI documentation with Doxygen"
            VERBATIM )
    else (DOXYGEN_FOUND)
        message("Doxygen need to be installed to generate WRF_HYDRO_NWM_JEDI documentation")
    endif (DOXYGEN_FOUND)
endif (BUILD_WRF_HYDRO_NWM_JEDI_DOC)
