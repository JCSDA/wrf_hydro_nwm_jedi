/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_TRAITS_H_
#define WRF_HYDRO_NWM_JEDI_TRAITS_H_

#include <string>

// #include "wrf_hydro_nwm_jedi/Covariance/Covariance.h"
#include "wrf_hydro_nwm_jedi/Geometry/Geometry.h"
// #include "wrf_hydro_nwm_jedi/GeometryIterator/GeometryIterator.h"
// #include "wrf_hydro_nwm_jedi/GetValuesTraj/GetValuesTraj.h"
// #include "wrf_hydro_nwm_jedi/Increment/Increment.h"
// #include "wrf_hydro_nwm_jedi/ModelAux/ModelAuxCovariance.h"
// #include "wrf_hydro_nwm_jedi/ModelAux/ModelAuxControl.h"
// #include "wrf_hydro_nwm_jedi/ModelAux/ModelAuxIncrement.h"
// #include "wrf_hydro_nwm_jedi/State/State.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasCovariance.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperator.h"

namespace wrf_hydro_nwm_jedi {

  struct Traits{
    static std::string name() {return "wrf_hydro_nwm_jedi";}
    static std::string nameCovar() {return "wrf_hydro_nwm_jediCovar";}
    static std::string nameCovar4D() {return "wrf_hydro_nwm_jediCovar";}

    // Interfaces that wrf_hydro_nwm_jedi has to implement
    // ---------------------------------------------------
    // typedef wrf_hydro_nwm_jedi::Covariance          Covariance;
    typedef wrf_hydro_nwm_jedi::Geometry            Geometry;
    // typedef wrf_hydro_nwm_jedi::GeometryIterator    GeometryIterator;
    // typedef wrf_hydro_nwm_jedi::GetValuesTraj       InterpolatorTraj;
    // typedef wrf_hydro_nwm_jedi::Increment           Increment;
    // typedef wrf_hydro_nwm_jedi::State               State;
    // typedef wrf_hydro_nwm_jedi::ModelAuxCovariance  ModelAuxCovariance;
    // typedef wrf_hydro_nwm_jedi::ModelAuxControl     ModelAuxControl;
    // typedef wrf_hydro_nwm_jedi::ModelAuxIncrement   ModelAuxIncrement;

    // Interfaces that are already provided by JEDI
    typedef ufo::GeoVaLs              GeoVaLs;
    typedef ufo::Locations            Locations;
    typedef ufo::ObsBias              ObsAuxControl;
    typedef ufo::ObsBiasCovariance    ObsAuxCovariance;
    typedef ufo::ObsBiasIncrement     ObsAuxIncrement;
    typedef ufo::ObsDiagnostics       ObsDiagnostics;
    typedef ufo::ObsOperator          ObsOperator;
    typedef ioda::ObsSpace            ObsSpace;
    typedef ioda::ObsVector           ObsVector;
    template <typename DATA> using ObsDataVector = ioda::ObsDataVector<DATA>;
  };
}  // namespace wrf_hydro_nwm_jedi

#endif  // WRF_HYDRO_NWM_JEDI_TRAITS_H_
