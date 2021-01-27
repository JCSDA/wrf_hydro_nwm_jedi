/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "wrf_hydro_nwm_jedi/Utilities/interface.h"

namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace wrf_hydro_nwm_jedi {

extern "C" {

  void wrf_hydro_nwm_jedi_getvalues_create_f90(F90getvalues &, const F90geom &, const F90locs &);

  void wrf_hydro_nwm_jedi_getvalues_delete_f90(F90getvalues &);

  void wrf_hydro_nwm_jedi_getvalues_fill_geovals_f90(
        const F90getvalues &, const F90geom &, const F90state &,
        const util::DateTime **, const util::DateTime **,
        const F90locs &, const F90goms &);

};  // extern "C"

// -------------------------------------------------------------------------------------------------

}  // namespace wrf_hydro_nwm_jedi
