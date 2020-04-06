/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_STATE_STATEFORTRAN_H_
#define WRF_HYDRO_NWM_JEDI_STATE_STATEFORTRAN_H_

#include "eckit/config/Configuration.h"

namespace wrf_hydro_nwm_jedi {

  typedef int F90state;

  extern "C" {
    void wrf_hydro_nwm_jedi_state_create_f90(F90state &, const F90geom &, const oops::Variables &);
    void wrf_hydro_nwm_jedi_state_read_file_f90(const F90geom &, const F90state &, const eckit::Configuration * const *, util::DateTime * const *);
    void wrf_hydro_nwm_jedi_state_copy_f90(const F90state &, const F90state &);
    /* void wrf_hydro_nwm_jedi_geometry_clone_f90(F90geom &, const F90geom &); */
    /* void wrf_hydro_nwm_jedi_geometry_delete_f90(F90geom &); */
  }
}  // namespace wrf_hydro_nwm_jedi

#endif  // WRF_HYDRO_NWM_JEDI_STATE_STATEFORTRAN_H_
