/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_UTILITIES_INTERFACE_H_
#define WRF_HYDRO_NWM_JEDI_UTILITIES_INTERFACE_H_

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace wrf_hydro_nwm_jedi {

// Geometry key type
typedef int F90geom;
// Geometry iterator key type
typedef int F90iter;
// Model key type
typedef int F90model;
// Tlm key type
typedef int F90tlm;
// Locations key type
typedef int F90locs;
// Goms key type
typedef int F90goms;
// Trajectory key type
typedef int F90traj;
// Background error covariance key type
typedef int F90bmat;
// Localization matrix
typedef int F90lclz;
// ObOp trajectory
typedef int F90ootrj;
// LinearVariableChange keys
typedef int F90lvcc2a;
typedef int F90lvca2m;
typedef int F90lvcnmcbal;
// VariableChange keys
typedef int F90vca2m;
typedef int F90vcc2a;
typedef int F90vcd2a;
typedef int F90m2g;
typedef int F90lm2g;
// State key
typedef int F90state;
// Increment key
typedef int F90inc;
// GetValues key
typedef int F90getvalues;
typedef int F90lineargetvalues;

}  // namespace wrf_hydro_nwm_jedi
#endif  // WRF_HYDRO_NWM_JEDI_UTILITIES_INTERFACE_H_
