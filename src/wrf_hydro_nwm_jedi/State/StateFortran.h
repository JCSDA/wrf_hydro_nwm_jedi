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
  typedef int F90inc;
  
  extern "C" {
    void wrf_hydro_nwm_jedi_state_create_f90(F90state &, const F90geom &, const oops::Variables &);
    void wrf_hydro_nwm_jedi_state_create_from_other_f90(F90state &, const F90state &);
    void wrf_hydro_nwm_jedi_state_read_file_f90(const F90geom &, const F90state &, const eckit::Configuration * const *, util::DateTime * const *);
    void wrf_hydro_nwm_jedi_state_print_f90(const F90state &, char *string);
    void wrf_hydro_nwm_jedi_state_get_mean_stddev_f90(const F90state &, int nf, float pstat[][1]);
    void wrf_hydro_nwm_jedi_state_copy_f90(const F90state &, const F90state &);
    void wrf_hydro_nwm_jedi_state_delete_f90(const F90state &);
    void wrf_hydro_nwm_jedi_state_add_incr_f90(const F90geom &, const F90state &, const F90inc &); 
    /* void wrf_hydro_nwm_jedi_geometry_clone_f90(F90geom &, const F90geom &); */
    /* void wrf_hydro_nwm_jedi_geometry_delete_f90(F90geom &); */
  }
}  // namespace wrf_hydro_nwm_jedi

#endif  // WRF_HYDRO_NWM_JEDI_STATE_STATEFORTRAN_H_
