/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_TRAITS_H_
#define WRF_HYDRO_NWM_JEDI_TRAITS_H_

#include <string>

#include "wrf_hydro_nwm_jedi/Covariance/Covariance.h"
#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
#include "wrf_hydro_nwm_jedi/GeometryIterator/GeometryIterator.h"
#include "wrf_hydro_nwm_jedi/Increment/Increment.h"
#include "wrf_hydro_nwm_jedi/LinearVariableChange/LinearVariableChange.h"
#include "wrf_hydro_nwm_jedi/ModelAux/ModelAuxControl.h"
#include "wrf_hydro_nwm_jedi/ModelAux/ModelAuxCovariance.h"
#include "wrf_hydro_nwm_jedi/ModelAux/ModelAuxIncrement.h"
#include "wrf_hydro_nwm_jedi/State/State.h"
#include "wrf_hydro_nwm_jedi/VariableChange/VariableChange.h"

// #include "wrf_hydro_nwm_jedi/Fields/Fields.h"


namespace wrf_hydro_nwm_jedi {

  struct Traits{
    static std::string name() {return "wrf_hydro_nwm_jedi";}
    static std::string nameCovar() {return "wrf_hydro_nwm_jediCovar";}
    static std::string nameCovar4D() {return "wrf_hydro_nwm_jediCovar";}

    // Interfaces that wrf_hydro_nwm_jedi has to implement
    // ---------------------------------------------------
    typedef wrf_hydro_nwm_jedi::Covariance          Covariance;
    typedef wrf_hydro_nwm_jedi::Geometry            Geometry;
    typedef wrf_hydro_nwm_jedi::GeometryIterator    GeometryIterator;
    typedef wrf_hydro_nwm_jedi::Increment             Increment;
    typedef wrf_hydro_nwm_jedi::State                 State;
    typedef wrf_hydro_nwm_jedi::ModelAuxCovariance    ModelAuxCovariance;
    typedef wrf_hydro_nwm_jedi::ModelAuxControl       ModelAuxControl;
    typedef wrf_hydro_nwm_jedi::ModelAuxIncrement     ModelAuxIncrement;
    typedef wrf_hydro_nwm_jedi::LinearVariableChange  LinearVariableChange;
    typedef wrf_hydro_nwm_jedi::VariableChange        VariableChange;
  };
}  // namespace wrf_hydro_nwm_jedi

#endif  // WRF_HYDRO_NWM_JEDI_TRAITS_H_
