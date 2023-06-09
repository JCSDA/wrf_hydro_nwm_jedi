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

  # GeometryIterator
  GeometryIterator/GeometryIterator.h
  GeometryIterator/GeometryIterator.cc
  GeometryIterator/GeometryIteratorFortran.h
  GeometryIterator/GeometryIterator.interface.F90
  GeometryIterator/GeometryIterator.F90
   
  # VariableChanges
  VariableChange/Model2GeoVaLs/VarChaModel2GeoVaLs.h
  VariableChange/Model2GeoVaLs/VarChaModel2GeoVaLs.cc
  VariableChange/Model2GeoVaLs/VarChaModel2GeoVaLsFortran.h
  VariableChange/Model2GeoVaLs/VarChaModel2GeoVaLs.interface.F90
  VariableChange/Model2GeoVaLs/VarChaModel2GeoVaLs.F90
  VariableChange/Base/VariableChangeBase.cc
  VariableChange/Base/VariableChangeBase.h
  VariableChange/VariableChange.cc
  VariableChange/VariableChange.h

  # LinearVariableChange
  LinearVariableChange/Model2GeoVaLs/LinVarChaModel2GeoVaLs.cc
  LinearVariableChange/Model2GeoVaLs/LinVarChaModel2GeoVaLs.h
  LinearVariableChange/Model2GeoVaLs/LinVarChaModel2GeoVaLs.F90
  LinearVariableChange/Model2GeoVaLs/LinVarChaModel2GeoVaLs.interface.F90
  LinearVariableChange/Model2GeoVaLs/LinVarChaModel2GeoVaLsFortran.h
  LinearVariableChange/Base/LinearVariableChangeBase.cc
  LinearVariableChange/Base/LinearVariableChangeBase.h
  LinearVariableChange/LinearVariableChange.cc
  LinearVariableChange/LinearVariableChange.h

  # # LinearGetValues
  #LinearGetValues/LinearGetValues.h
  #LinearGetValues/LinearGetValues.cc
   
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

  # ObsLocalization
  ObsLocalization/instantiateObsLocFactory.h
  ObsLocalization/ObsLocVerticalBrasnett.h

  # Utilities
  Utilities/utilities.F90
  Utilities/constants.F90	
  )


ecbuild_add_library( TARGET   wrf_hydro_nwm_jedi
                     SOURCES  ${wrf_hydro_nwm_jedi_src_files}
                     LIBS     ${NETCDF_LIBRARIES}
                              eckit eckit_mpi fckit oops saber ioda ufo
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${WRF_HYDRO_NWM_JEDI_LINKER_LANGUAGE}
                    )
