#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

list( APPEND wrf_hydro_nwm_jedi_src_files
  Traits.h
  
  # # Covariance
  Covariance/Covariance.h
  Covariance/Covariance.cc
  Covariance/CovarianceFortran.h
  Covariance/Covariance.interface.F90
  Covariance/Covariance.F90
  
  # Geometry  
  Geometry/Geometry.h  
  Geometry/Geometry.cc
  Geometry/GeometryFortran.h
  Geometry/Geometry.interface.F90
  Geometry/Geometry.F90

  # # Fields
  Fields/Fields.h
  Fields/Fields.cc
  Fields/Fields.F90

  # # GeometryIterator
  # GeometryIterator/GeometryIterator.h
  # GeometryIterator/GeometryIterator.cc

  # GetValues
  GetValues/GetValues.h
  GetValues/GetValues.cc
  GetValues/GetValues.F90
  GetValues/LinearGetValues.cc
  GetValues/LinearGetValues.h
  GetValues/GetValues.interface.F90
  GetValues/LinearGetValues.interface.h
  GetValues/LinearGetValues.interface.F90
  
  # VariableChanges
  VariableChanges/Model2GeoVaLs/VarChaModel2GeoVaLs.cc
  VariableChanges/Model2GeoVaLs/VarChaModel2GeoVaLs.h

  # # LinearGetValues
  #LinearGetValues/LinearGetValues.h
  #LinearGetValues/LinearGetValues.cc

  # # GetValuesTraj
  # GetValuesTraj/GetValuesTraj.h
  # GetValuesTraj/GetValuesTraj.cc
  
  # # Increment
  Increment/Increment.interface.F90
  Increment/Increment.h
  Increment/Increment.cc
  Increment/Increment_utils.F90
  Increment/IncrementFortran.h
  Increment/Increment.F90

  # LinearModel
  LinearModel/TlmId.h
  LinearModel/TlmId.cc

  # Model
  Model/Model.h
  Model/Model.cc

  # ModelAux
  ModelAux/ModelAuxControl.h
  ModelAux/ModelAuxIncrement.h
  ModelAux/ModelAuxCovariance.h
  
  # State
  State/State.h
  State/State.cc
  State/State_utils.F90
  State/State.interface.F90
  State/StateFortran.h
  State/State.F90

  # Utilities
  Utilities/utilities.F90
  Utilities/constants.F90	
  )


ecbuild_add_library( TARGET   wrf_hydro_nwm_jedi
                     SOURCES  ${wrf_hydro_nwm_jedi_src_files}
                     LIBS     ${NETCDF_LIBRARIES}
                              eckit eckit_mpi fckit oops ioda ufo
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${WRF_HYDRO_NWM_JEDI_LINKER_LANGUAGE}
                    )
