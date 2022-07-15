/*
 * (C) Copyright 2019-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef WRF_HYDRO_NWM_JEDI_STATE_STATEFORTRAN_H_
#define WRF_HYDRO_NWM_JEDI_STATE_STATEFORTRAN_H_

#include <string>

#include "eckit/config/Configuration.h"

namespace wrf_hydro_nwm_jedi {

  typedef int F90state;
  typedef int F90inc;

  extern "C" {
    void wrf_hydro_nwm_jedi_state_create_f90(
        F90state &, const F90geom &, const oops::Variables &);
    void wrf_hydro_nwm_jedi_state_create_from_other_f90(
        F90state &, const F90state &);
    void wrf_hydro_nwm_jedi_state_read_file_f90(
        const F90geom &, const F90state &,
        const eckit::Configuration * const *,
        const util::DateTime * const *);
    void wrf_hydro_nwm_jedi_state_write_file_f90(
        const F90geom &,
        const F90state &,
        const eckit::Configuration * const *,
        const util::DateTime * const *);
    void wrf_hydro_nwm_jedi_state_print_f90(const F90state &, char *string);
    // void wrf_hydro_nwm_jedi_state_get_mean_stddev_f90(
    //     const F90state &, int nf, float pstat[][1]);
    void wrf_hydro_nwm_jedi_state_copy_f90(const F90state &, const F90state &);
    void wrf_hydro_nwm_jedi_state_delete_f90(const F90state &);
    void wrf_hydro_nwm_jedi_state_add_incr_f90(const F90state &, const F90inc &);
    double wrf_hydro_nwm_jedi_state_rms_f90(const F90state &);
    void wrf_hydro_nwm_jedi_state_zero_f90(const F90inc &);
    void wrf_hydro_nwm_jedi_state_ones_f90(const F90inc &);
    void wrf_hydro_nwm_jedi_state_axpy_f90(
        const F90state &, const double &, const F90state &);
    void wrf_hydro_nwm_jedi_state_to_fieldset_f90(const F90state &, const F90geom &,
                  const oops::Variables &, atlas::field::FieldSetImpl *);
    void wrf_hydro_nwm_jedi_state_from_fieldset_f90(const F90state &, const F90geom &,
                  const oops::Variables &, const atlas::field::FieldSetImpl *);
  }

}  // namespace wrf_hydro_nwm_jedi

#endif  // WRF_HYDRO_NWM_JEDI_STATE_STATEFORTRAN_H_
