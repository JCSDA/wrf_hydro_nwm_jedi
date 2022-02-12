/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "wrf_hydro_nwm_jedi/Utilities/interface.h"

namespace wrf_hydro_nwm_jedi {
  typedef int F90lvc_M2G;
  extern "C" {
  void wrf_hydro_nwm_jedi_lvc_model2geovals_create_f90(const F90lvc_M2G &, const F90geom &,
                                            const F90state &, const F90state &,
                                            const eckit::LocalConfiguration &);
  void wrf_hydro_nwm_jedi_lvc_model2geovals_delete_f90(F90lvc_M2G &);
  void wrf_hydro_nwm_jedi_lvc_model2geovals_multiply_f90(const F90lvc_M2G &, const F90geom &,
                                            const F90inc &, const F90inc &);
  void wrf_hydro_nwm_jedi_lvc_model2geovals_multiplyadjoint_f90(const F90lvc_M2G &,
                                            const F90geom &, const F90inc &, const F90inc &);
  }  // extern "C"
}  // namespace wrf_hydro_nwm_jedi
