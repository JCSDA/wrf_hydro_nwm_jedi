/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "wrf_hydro_nwm_jedi/Utilities/interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace wrf_hydro_nwm_jedi {
  typedef int F90vc_M2G;
  extern "C" {
  void wrf_hydro_nwm_jedi_varchamodel2geovals_create_f90(const F90vc_M2G &, const F90geom &,
                                           const eckit::Configuration * const *);
  void wrf_hydro_nwm_jedi_varchamodel2geovals_delete_f90(F90vc_M2G &);
  void wrf_hydro_nwm_jedi_varchamodel2geovals_changevar_f90(const F90vc_M2G &, const F90geom &,
                                               const F90state &, const F90state &);
}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm_jedi
